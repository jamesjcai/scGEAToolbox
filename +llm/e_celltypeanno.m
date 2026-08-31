function [c_cell_type_tx, cL_cell_type_tx, T] = e_celltypeanno(sce, options)
%E_CELLTYPEANNO Annotate clusters by asking a local LLM about their markers.
%
%   [c_cell_type_tx, cL_cell_type_tx, T] = LLM.E_CELLTYPEANNO(sce)
%   [...] = LLM.E_CELLTYPEANNO(sce, Name=Value)
%
%   Ranks each cluster's marker genes (PKG.E_FINDALLMARKERS), asks a locally
%   hosted model to name the cell type they describe, and returns the labels
%   in the same three-output shape as SC_CELLTYPEANNO, so the two are
%   interchangeable at the call site.
%
%   The two differ in kind, not quality: SC_CELLTYPEANNO matches markers
%   against PanglaoDB through RUN.ML_ALONA_NEW and is deterministic; this
%   asks a language model, which is not. See REPRODUCIBILITY below.
%
%   OUTPUTS
%     c_cell_type_tx  - NumCells-by-1 string, one label per cell
%     cL_cell_type_tx - one label per cluster, in cluster order
%     T               - table: Cluster, CellType, Agreement, NumVotes,
%                       Markers (the genes the model was shown)
%
%   NAME-VALUE ARGUMENTS
%     Species      "human" (default), "mouse", or any species name
%     Tissue       tissue or organ, e.g. "intestine". Default "" (unstated).
%                  Worth setting - marker sets are far less ambiguous once
%                  the tissue is known.
%     NumGenes     marker genes shown per cluster. Default 30.
%     Provider     LLM provider. Default: whatever gui.i_setllmmodel
%                  configured. Any provider with a client in +llm works -
%                  Ollama, OpenAI, TAMUAIChat, NVIDIA, Gemini - so this does
%                  not require a local model.
%     Model        model name for that provider. Default: the configured one.
%     NIterations  independent answers per cluster, majority-voted.
%                  Default 3. See REPRODUCIBILITY.
%     ClusterID    grouping vector to use instead of sce.c_cluster_id
%     Markers      a precomputed PKG.E_FINDALLMARKERS table, to skip the DE step
%     Verbose      print progress. Default true.
%
%   REPRODUCIBILITY
%   A language model is sampled, not evaluated: the same cluster can get
%   different names on different runs, and nothing here makes that go away.
%   NIterations asks independently several times and takes the majority,
%   which stabilises the answer and - more usefully - measures its own
%   instability. T.Agreement is the winning label's share of the votes, so a
%   cluster answered "T cells" three times out of three (1.00) and one
%   answered three different ways (0.33) are distinguishable. Treat a low
%   Agreement as the model telling you it does not know.
%
%   The model is also told to answer "Unknown" rather than guess. Whether it
%   obeys is a property of the model, not of this function.
%
%   REQUIREMENTS
%   A configured LLM provider - set one with GUI.I_SETLLMMODEL (menu Options
%   > Set LLM Provider & Model...). A hosted provider needs only its API key
%   in the configured .env; a local one needs Ollama running, for which
%   LLM.E_INSTALLOLLAMA can install a portable copy.
%
%   The provider is smoke-tested with one short prompt BEFORE the marker
%   step, so a missing key, an unreachable server or an uninstalled model
%   surfaces as a single clear error rather than one failure per cluster
%   after the expensive part has already run.
%
%   EXAMPLE
%     sce = sce.clustercells(12);
%     [ct, cL, T] = llm.e_celltypeanno(sce, Species="mouse", Tissue="pancreas");
%     sce.c_cell_type_tx = ct;
%     disp(T(:, ["Cluster" "CellType" "Agreement"]));
%
%   See also SC_CELLTYPEANNO, PKG.E_FINDALLMARKERS, GUI.I_SETLLMMODEL,
%   LLM.I_ASKLLM, LLM.I_LLMSETTINGS.

arguments
    sce (1,1) SingleCellExperiment
    options.Species (1,1) string = "human"
    options.Tissue (1,1) string = ""
    options.NumGenes (1,1) double {mustBeInteger, mustBePositive} = 30
    options.Model (1,1) string = ""
    options.Provider (1,1) string = ""
    options.NIterations (1,1) double {mustBeInteger, mustBePositive} = 3
    options.ClusterID = []
    options.Markers = []
    options.Verbose (1,1) logical = true
end

% Every failure below raises rather than returning early, so the outputs are
% assigned on every path that reaches the end.

% ---- clusters -----------------------------------------------------------
if isempty(options.ClusterID)
    rawc = sce.c_cluster_id;
else
    rawc = options.ClusterID;
end
if isempty(rawc)
    error('llm:e_celltypeanno:NoClusters', ...
        'No clustering found. Run sce.clustercells first, or pass ClusterID.');
end
if numel(rawc) ~= sce.NumCells
    error('llm:e_celltypeanno:ClusterLength', ...
        'ClusterID has %d entries for %d cells.', numel(rawc), sce.NumCells);
end
% Natural order, so cluster 10 comes after cluster 9 rather than after
% cluster 1. findgroups(string(...)) alone sorts lexicographically.
[cidx, cL] = pkg.i_grp2idxsorted(rawc(:));
nClust = numel(cL);

% ---- preflight ----------------------------------------------------------
% Check before computing markers: the DE step is the expensive part, and
% failing after it because no model is reachable wastes all of it.
[provider, model] = llm.i_llmsettings(options.Model, options.Provider);

% The reference dictionary is per-tissue, so it is fetched once here rather
% than rebuilt for every cluster. It comes back empty when no tissue was
% given, when the tissue is not catalogued, or when the download failed -
% all of which just mean the prompt goes out without a reference section.
[ctx, ctxtypes] = llm.i_cellcontext(options.Tissue);
if options.Verbose && strlength(ctx) > 0
    fprintf('[celltypeanno] reference dictionary: %d cell type(s) for %s\n', ...
        numel(ctxtypes), options.Tissue);
end
i_smoketest(provider, model);
if options.Verbose
    fprintf('[celltypeanno] %d cluster(s), %s / %s, %d vote(s) each\n', ...
        nClust, provider, model, options.NIterations);
end

% ---- markers ------------------------------------------------------------
if isempty(options.Markers)
    if options.Verbose, fprintf('[celltypeanno] ranking markers...\n'); end
    % A generous cap: e_findallmarkers truncates each cluster to
    % maxnummarkers BEFORE sorting by significance, so a small cap would
    % select genes by their position in sce.g rather than by how well they
    % mark the cluster. Ask for many and rank them here instead.
    Tm = pkg.e_findallmarkers(sce.X, sce.g, cidx, cL, 0.5, 0.1, false, 500);
else
    Tm = options.Markers;
end
if isempty(Tm) || height(Tm) == 0
    error('llm:e_celltypeanno:NoMarkers', ...
        'No marker genes passed the significance filter for any cluster.');
end

% ---- ask ----------------------------------------------------------------
labels = strings(nClust, 1);
agree = zeros(nClust, 1);
votes = zeros(nClust, 1);
shown = strings(nClust, 1);

for k = 1:nClust
    genes = i_markersfor(Tm, cL(k), options.NumGenes);
    shown(k) = strjoin(genes, ", ");
    if isempty(genes)
        labels(k) = "Unknown";
        continue
    end

    prompt = i_prompt(options.Species, options.Tissue, cL(k), genes, ctx);
    answers = strings(0, 1);
    for it = 1:options.NIterations
        [raw, ok] = llm.i_askllm(prompt, provider, model);
        if ~ok
            warning('llm:e_celltypeanno:GenerateFailed', ...
                'Cluster %s, attempt %d: %s', cL(k), it, raw);
            continue
        end
        lbl = llm.i_cleanllmlabel(raw);
        if strlength(lbl) > 0
            answers(end+1) = lbl; %#ok<AGROW>
        end
    end

    [labels(k), agree(k), votes(k)] = llm.i_votelabels(answers);
    if options.Verbose
        fprintf('[celltypeanno] cluster %-12s -> %-28s (agreement %.2f, %d vote(s))\n', ...
            cL(k), labels(k), agree(k), votes(k));
    end
end

% ---- assemble -----------------------------------------------------------
cL_cell_type_tx = labels;
c_cell_type_tx = labels(cidx);
T = table(cL(:), labels, agree, votes, shown, 'VariableNames', ...
    {'Cluster', 'CellType', 'Agreement', 'NumVotes', 'Markers'});
end


% =========================================================================
function i_smoketest(provider, model)
% One short prompt before the expensive marker step. A missing API key, an
% unreachable server or an uninstalled model all surface here, once, instead
% of once per cluster per vote after the DE has already been computed.
[txt, ok] = llm.i_askllm("Reply with the single word: OK", provider, model);
if ~ok
    error('llm:e_celltypeanno:ProviderUnavailable', ...
        ['The configured LLM (%s / %s) did not answer a test prompt:\n  %s\n' ...
         'Check the provider and model with gui.i_setllmmodel, and that the ' ...
         'API key file it points at is valid.'], provider, model, txt);
end
end


function genes = i_markersfor(Tm, cluster, n)
% The top n markers for one cluster. e_findallmarkers returns rows sorted by
% p_val_adj within each group, so taking the head is taking the strongest.
rows = string(Tm.grp) == string(cluster);
if ~any(rows)
    genes = strings(0, 1);
    return
end
sub = Tm(rows, :);
if ismember('p_val_adj', sub.Properties.VariableNames)
    sub = sortrows(sub, 'p_val_adj', 'ascend');
end
genes = string(sub.g(1:min(n, height(sub))));
genes = genes(:)';
end


function p = i_prompt(species, tissue, cluster, genes, ctx)
context = "";
if strlength(tissue) > 0
    context = " from " + tissue;
end

% With a reference dictionary the task changes: the model picks from a
% fixed vocabulary rather than naming a cell type from memory. Say so
% explicitly, and allow Unknown so it is not forced into a bad match.
reference = "";
if strlength(ctx) > 0
    reference = "CELL MARKERS FOR CONTEXT (Reference Dictionary):" + newline + ...
        ctx + newline + newline + ...
        "Prefer a cell type named in the reference list above, choosing the " + ...
        "most specific one whose markers overlap the cluster's markers. " + ...
        "Only name a cell type outside the list if none of them fit." + ...
        newline + newline;
end

p = "You are annotating single-cell RNA-seq clusters." + newline + newline + ...
    "Species: " + species + newline + ...
    "Tissue: " + i_orunstated(tissue) + newline + ...
    "Cluster: " + cluster + newline + newline + ...
    reference + ...
    "Top marker genes for this cluster, most significant first:" + newline + ...
    strjoin(genes, ", ") + newline + newline + ...
    "Name the single most likely cell type for this cluster of cells" + ...
    context + ". Reply with the cell type name ONLY - no explanation, " + ...
    "no preamble, no punctuation, no markdown. If the markers do not " + ...
    "identify a cell type confidently, reply exactly: Unknown";
end


function s = i_orunstated(tissue)
if strlength(tissue) > 0
    s = tissue;
else
    s = "not stated";
end
end
