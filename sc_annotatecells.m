function [sce, T, stashname, methodResult] = sc_annotatecells(sce, options)
%SC_ANNOTATECELLS Assign cell types by any available method, through one call.
%
%   [sce, T, stashname] = SC_ANNOTATECELLS(sce)
%   [sce, T] = SC_ANNOTATECELLS(sce, Method="llm", Tissue="pancreas")
%   [sce, T] = SC_ANNOTATECELLS(sce, Method="consensus", Methods=["markers" "llm"])
%
%   Routes to one of the toolbox's cell-type annotation methods, writes the
%   result to sce.c_cell_type_tx, and returns a table describing what was
%   assigned. The individual entry points keep working unchanged; this is a
%   front door, not a replacement.
%
%   METHODS
%     "markers"     (default) marker matching against PanglaoDB, via
%                   SC_CELLTYPEANNO. Deterministic. Needs clusters.
%                   The FDR column gives a hypergeometric test of the
%                   chosen type's markers against the cluster's own.
%     "llm"         a language model reading each cluster's marker genes, via
%                   LLM.E_CELLTYPEANNO. Needs clusters and a configured
%                   provider (GUI.I_SETLLMMODEL). Not deterministic - see
%                   Confidence below.
%     "scimilarity" reference-model annotation, via RUN.PY_SCIMILARITY.
%                   Per-cell, no clustering needed. Needs ModelDir and a
%                   working Python environment.
%     "panhumanpy"  reference-model annotation, via RUN.PY_PANHUMANPY.
%                   Per-cell, no clustering needed. Needs Python.
%     "custommarkers"  marker matching against a list you supply in
%                   Markers. Each marker is scored against the other
%                   clusters, so an abundant transcript present in every
%                   cell (ambient insulin, say) cannot win everywhere.
%                   Needs clusters. A cluster expressing no list's markers
%                   above the other clusters is "undetermined". FDR as for
%                   "markers".
%     "reference"   comparison with a summary-statistic reference - your
%                   own atlas or a published one - via SC_CELLTYPEANNOREF.
%                   Needs clusters and Reference (see SC_BUILDCELLTYPEREF).
%                   No raw reference data and no integration.
%     "consensus"   run each of Methods and keep the label they agree on.
%                   Needs clusters, because agreement is judged per cluster.
%
%   THE METHODS ARE NOT INTERCHANGEABLE
%   A common interface makes them easy to swap, which makes it easy to forget
%   they answer differently. "markers" and "llm" label CLUSTERS, so they
%   inherit whatever the clustering got wrong. "scimilarity" and "panhumanpy"
%   label CELLS against a reference, so they need no clustering but do need a
%   multi-gigabyte model. They also differ in granularity: on a pancreas
%   dataset whose ground truth said "Enteroendocrine cells", the LLM answered
%   "Delta cells" - finer, and not wrong. T.Method is therefore part of the
%   result, not bookkeeping: a label means little without knowing what
%   produced it.
%
%   OUTPUT TABLE
%     Method     which method produced the row
%     Cluster    the cluster labelled, or "" for the per-cell methods
%     CellType   the assigned label
%     NumCells   how many cells carry it
%     Confidence 0-1 where the method reports one, NaN where it does not.
%                Only "llm" (vote agreement) and "consensus" (share of
%                methods agreeing) can fill this. It is left NaN rather than
%                set to 1, because "deterministic" is not "certain" and a
%                fabricated confidence is worse than an absent one.
%     FDR        BH-adjusted hypergeometric p-value that the cluster's
%                markers overlap the assigned type's markers, for
%                "markers" and "reference"; NaN otherwise, or where the
%                cluster had no marker passing the filter of
%                PKG.I_CLUSTERMARKERS. A small FDR supports the label; a
%                large one says the markers do not, whatever produced it.
%
%   With two methods, a Confidence of 0.5 means they disagreed outright and
%   NO majority exists. The label reported is then the earlier entry in
%   Methods - a documented precedence, not a judgement about which is right.
%   Filter on Confidence < 1 to see every cluster the methods split on.
%
%   NAME-VALUE ARGUMENTS
%     Method, Methods, Species, Tissue, ModelDir, WorkDir,
%     Provider, Model, NIterations  - passed to the method that uses them
%     Reference        reference struct for "reference"
%     Markers          marker list for "custommarkers": a two-column table,
%                      cell type name then its genes as one comma separated
%                      string - what GUI.I_GETCUSTOMMARKERS and
%                      PKG.I_PARSEMARKERLIST return
%     ReferenceMethod  "both" (default), "correlation" or "markers";
%                      passed to SC_CELLTYPEANNOREF as Method
%     KeepOld    stash any existing labels in the 'old_cell_type' cell
%                attribute before overwriting. Default true.
%     Verbose    print progress. Default true.
%
%   STASHNAME is the 'old_cell_type_N' cell attribute that the labels being
%   replaced were kept in, or "" if there were none. A GUI caller needs it
%   to tell the user where they went; see GUI.I_STASHNOTICE.
%
%   METHODRESULT is the method's own, fuller result where it has one that
%   the common table leaves out - for "reference", the table from
%   SC_CELLTYPEANNOREF, with both tests' scores and whether they agree -
%   and [] otherwise.
%
%   PROVENANCE
%   The method, its settings and the time are appended to sce.metadata, so a
%   saved object can still say whether its labels came from PanglaoDB or from
%   a language model. Nothing else in the toolbox records that today.
%
%   EXAMPLE
%     sce = sce.clustercells(12);
%     [sce, T] = sc_annotatecells(sce, Method="consensus", ...
%                                 Methods=["markers" "llm"], Species="mouse");
%     disp(T(T.Confidence < 1, :));      % where the methods disagreed
%
%   See also SC_CELLTYPEANNO, SC_CELLTYPEANNOREF, LLM.E_CELLTYPEANNO, RUN.PY_SCIMILARITY,
%   RUN.PY_PANHUMANPY, SINGLECELLEXPERIMENT/ASSIGNCELLTYPE.

arguments
    sce (1,1) SingleCellExperiment
    options.Method (1,1) string {mustBeMember(options.Method, ...
        ["markers", "llm", "scimilarity", "panhumanpy", "reference", ...
        "custommarkers", "consensus"])} = "markers"
    options.Methods (1,:) string {mustBeMember(options.Methods, ...
        ["markers", "llm", "scimilarity", "panhumanpy", "reference", ...
        "custommarkers"])} = ...
        ["markers", "llm"]
    options.Species (1,1) string = "human"
    options.Tissue (1,1) string = ""
    options.ModelDir (1,1) string = ""
    options.WorkDir (1,1) string = ""
    options.Provider (1,1) string = ""
    options.Model (1,1) string = ""
    options.NIterations (1,1) double {mustBeInteger, mustBePositive} = 3
    options.Reference struct = struct([])
    options.Markers table = table()
    options.ReferenceMethod (1,1) string {mustBeMember(options.ReferenceMethod, ...
        ["both", "correlation", "markers"])} = "both"
    options.KeepOld (1,1) logical = true
    options.Verbose (1,1) logical = true
end

stashname = "";
methodResult = [];

if options.Method == "consensus"
    [labels, T] = i_consensus(sce, options);
    detail = "methods=" + strjoin(options.Methods, "+");
else
    [labels, T, methodResult] = i_runone(sce, options.Method, options);
    detail = i_detailof(options.Method, options);
end

if isempty(labels)
    if options.Verbose
        fprintf('[annotate] no labels produced; sce is unchanged.\n');
    end
    return
end

% Overwriting labels silently loses them. The GUI annotation callbacks
% stash the previous ones the same way; this makes that the rule rather
% than the majority behaviour. Each call adds a new numbered attribute
% rather than replacing one slot, so a history of past annotations
% survives, not just the most recent one. Callbacks that repurpose
% c_cell_type_tx as a working group vector (GUI.CALLBACK_SCTENIFOLDCKO,
% GUI.CALLBACK_SCTENIFOLDXCT) deliberately do not stash.
if options.KeepOld
    stashname = pkg.i_stashcelltypehistory(sce);
    if stashname ~= "" && options.Verbose
        fprintf('[annotate] previous labels kept in the ''%s'' attribute.\n', ...
            stashname);
    end
end

sce.c_cell_type_tx = string(labels(:));
sce = i_recordprovenance(sce, options.Method, detail);

if options.Verbose
    fprintf('[annotate] %s: %d cell type(s) over %d cells\n', ...
        options.Method, numel(unique(string(labels))), sce.NumCells);
end
end


% =========================================================================
function [labels, T, methodResult] = i_runone(sce, method, options)
% One method, normalised to per-cell labels plus the common table.
%
% There are no LABELS/T defaults before the switch any more: every case
% assigns both and the OTHERWISE errors, so the defaults were only ever
% reachable by an unrecognised method name -- which is exactly the silent
% empty result that made a broken consensus look unanimous. METHODRESULT
% is optional extra detail, so [] is a real answer and it keeps a default.
methodResult = [];
switch method
    case "markers"
        i_requireclusters(sce, method);
        [c, cL] = pkg.i_grp2idxsorted(sce.c_cluster_id);
        [~, cLtx] = sc_celltypeanno(sce.X, sce.g, c, char(options.Species));
        cLtx = string(cLtx(:));
        labels = cLtx(c);
        [types, lists] = i_panglaodblists(options.Species);
        fdr = i_listfdr(sce, c, cLtx, types, lists);
        T = i_clustertable(method, cL, cLtx, c, nan(numel(cL), 1), fdr);

    case "custommarkers"
        i_requireclusters(sce, method);
        if isempty(options.Markers) || width(options.Markers) < 2
            error('sc_annotatecells:NoMarkers', ...
                ['Method "custommarkers" needs Markers, a two-column table ', ...
                 'of cell type names and comma separated marker genes.']);
        end
        [c, cL] = pkg.i_grp2idxsorted(sce.c_cluster_id);
        [cLtx, types, lists] = i_scorecustom(sce, c, numel(cL), options.Markers);
        labels = cLtx(c);
        fdr = i_listfdr(sce, c, cLtx, types, lists);
        T = i_clustertable(method, cL, cLtx, c, nan(numel(cL), 1), fdr);

    case "llm"
        i_requireclusters(sce, method);
        [c, cL] = pkg.i_grp2idxsorted(sce.c_cluster_id);
        [labels, cLtx, Tllm] = llm.e_celltypeanno(sce, ...
            Species = options.Species, Tissue = options.Tissue, ...
            Provider = options.Provider, Model = options.Model, ...
            NIterations = options.NIterations, Verbose = options.Verbose);
        labels = string(labels(:));
        T = i_clustertable(method, cL, string(cLtx(:)), c, Tllm.Agreement);

    case "scimilarity"
        if strlength(options.ModelDir) == 0
            error('sc_annotatecells:NoModelDir', ...
                ['Method "scimilarity" needs ModelDir pointing at a ', ...
                 'downloaded SCimilarity model.']);
        end
        wk = i_workdir(options.WorkDir);
        labels = string(run.py_scimilarity(sce, char(options.ModelDir), wk, '', true));
        T = i_celltable(method, labels);

    case "panhumanpy"
        wk = i_workdir(options.WorkDir);
        labels = string(run.py_panhumanpy(sce, wk, true));
        T = i_celltable(method, labels);

    case "reference"
        i_requireclusters(sce, method);
        if isempty(options.Reference)
            error('sc_annotatecells:NoReference', ...
                ['Method "reference" needs Reference, a struct made by ', ...
                 'SC_BUILDCELLTYPEREF or assembled from published tables.']);
        end
        [c, cL] = pkg.i_grp2idxsorted(sce.c_cluster_id);
        [cLtx, Tref] = sc_celltypeannoref(sce.X, sce.g, sce.c_cluster_id, ...
            options.Reference, Method=options.ReferenceMethod);
        labels = cLtx(c);
        T = i_clustertable(method, cL, cLtx, c, nan(numel(cL), 1), Tref.FDR);
        methodResult = Tref;
    otherwise
        % LABELS used to be preset empty above and this switch had no
        % OTHERWISE, so an unrecognised name returned no labels, which
        % i_consensus read as a method that had declined to answer.
        error('sc_annotatecells:UnknownMethod', ...
            ['Unknown method "%s". Use "markers", "llm", "scimilarity", ', ...
            '"panhumanpy", "reference" or "custommarkers".'], method);
end

if ~isempty(labels) && numel(labels) ~= sce.NumCells
    error('sc_annotatecells:LabelCount', ...
        ['Method "%s" returned %d labels for %d cells. The annotation was ', ...
         'not applied.'], method, numel(labels), sce.NumCells);
end
end


function [labels, T] = i_consensus(sce, options)
% Agreement between methods, judged per cluster. A per-cell method's vote for
% a cluster is the label most of that cluster's cells received, so every
% method can be compared on the same footing.
i_requireclusters(sce, "consensus");
if numel(options.Methods) < 2
    error('sc_annotatecells:ConsensusNeedsTwo', ...
        'Method "consensus" needs at least two entries in Methods.');
end

[c, cL] = pkg.i_grp2idxsorted(sce.c_cluster_id);
nClust = numel(cL);
votes = strings(nClust, numel(options.Methods));

for m = 1:numel(options.Methods)
    if options.Verbose
        fprintf('[annotate] consensus %d/%d: %s\n', m, numel(options.Methods), ...
            options.Methods(m));
    end
    [lab, ~] = i_runone(sce, options.Methods(m), options);
    if isempty(lab)
        % Not silently. The per-cell backends return empty on failure
        % rather than raising -- run.py_panhumanpy and run.py_scimilarity
        % both preset celltypes = [] and only fill it when the Python side
        % exits cleanly -- so without this the caller cannot tell a
        % unanimous two-method consensus from one method that never ran.
        warning('sc_annotatecells:methodProducedNoLabels', ...
            ['Method "%s" returned no labels, so it casts no vote. The ', ...
            'Confidence column counts it as a method that did not ', ...
            'agree.'], options.Methods(m));
        votes(:, m) = "";
        continue
    end
    for k = 1:nClust
        votes(k, m) = i_majority(string(lab(c == k)));
    end
end

cLtx = strings(nClust, 1);
conf = zeros(nClust, 1);
for k = 1:nClust
    [cLtx(k), conf(k)] = i_agree(votes(k, :));
end
labels = cLtx(c);
T = i_clustertable("consensus", cL, cLtx, c, conf);
end


% =========================================================================
function i_requireclusters(sce, method)
if isempty(sce.c_cluster_id) || isscalar(unique(sce.c_cluster_id))
    error('sc_annotatecells:NoClusters', ...
        ['Method "%s" labels clusters, so it needs more than one. Run ', ...
         'sce.clustercells first, or use a per-cell method ', ...
         '("scimilarity", "panhumanpy").'], method);
end
end


function wk = i_workdir(requested)
if strlength(requested) > 0
    wk = char(requested);
    if ~isfolder(wk), mkdir(wk); end
else
    wk = pkg.i_tempdirfile();
end
end


function lbl = i_majority(v)
v = v(strlength(v) > 0);
if isempty(v)
    lbl = "";
    return
end
[u, ~, idx] = unique(v);
[~, w] = max(accumarray(idx, 1));
lbl = u(w);
end


function [best, conf] = i_agree(v)
% The label most methods gave, and the share that gave it. This is the same
% majority-with-agreement-score that LLM.I_VOTELABELS already performs over
% repeated model answers, so it delegates rather than keeping a second copy
% that could drift - votes arrive in Methods order, and that helper breaks
% ties in favour of the first, which is the precedence documented above.
%
% Its "no answers" sentinel is "Unknown"; the toolbox's word for an unlabelled
% cell is "undetermined", so that one value is translated.
%
% The denominator is the number of methods ASKED FOR, which is numel(v).
% LLM.I_VOTELABELS strips blank answers before dividing -- right for
% repeated model answers, where a blank reply is not a vote to count against
% anything -- but here a blank means one of the requested methods produced
% nothing at all. Dividing by the methods that answered made a
% single-answer consensus look unanimous: with Methods of two, one of them
% silently returning empty gave Confidence 1.00, where the help says 1.00
% is agreement and 0.50 is outright disagreement, and describes the column
% as the "share of methods agreeing".
nRequested = numel(v);
[best, conf, n] = llm.i_votelabels(v);
if n == 0
    best = "undetermined";
elseif n < nRequested
    % conf is topCount/n; rescale to topCount/nRequested.
    conf = conf*n/nRequested;
end
end



function T = i_clustertable(method, cL, cLtx, c, conf, fdr)
n = numel(cL);
if nargin < 6 || numel(fdr) ~= n
    fdr = nan(n, 1);
end
counts = zeros(n, 1);
for k = 1:n
    counts(k) = sum(c == k);
end
if numel(conf) ~= n
    conf = nan(n, 1);
end
T = table(repmat(string(method), n, 1), string(cL(:)), string(cLtx(:)), ...
    counts, double(conf(:)), double(fdr(:)), 'VariableNames', ...
    {'Method', 'Cluster', 'CellType', 'NumCells', 'Confidence', 'FDR'});
end


function T = i_celltable(method, labels)
% Per-cell methods have no cluster to report, so Cluster stays empty rather
% than being filled with something invented.
[u, ~, idx] = unique(string(labels));
counts = accumarray(idx, 1);
n = numel(u);
T = table(repmat(string(method), n, 1), repmat("", n, 1), u(:), ...
    counts(:), nan(n, 1), nan(n, 1), 'VariableNames', ...
    {'Method', 'Cluster', 'CellType', 'NumCells', 'Confidence', 'FDR'});
T = sortrows(T, 'NumCells', 'descend');
end


function d = i_detailof(method, options)
switch method
    case "markers"
        d = "database=panglaodb, species=" + options.Species;
    case "llm"
        [p, m] = llm.i_llmsettings(options.Model, options.Provider);
        d = "provider=" + p + ", model=" + m + ...
            ", votes=" + options.NIterations + ", species=" + options.Species;
        if strlength(options.Tissue) > 0
            d = d + ", tissue=" + options.Tissue;
        end
    case "scimilarity"
        d = "model=" + options.ModelDir;
    case "reference"
        d = "reference=" + i_refname(options.Reference) + ...
            ", test=" + options.ReferenceMethod;
    case "custommarkers"
        d = "custom markers, " + height(options.Markers) + " cell types";
    otherwise
        d = "species=" + options.Species;
end
end


function name = i_refname(ref)
name = "custom";
if isfield(ref, 'Source') && strlength(string(ref.Source)) > 0
    name = string(ref.Source);
end
end


function [types, lists] = i_panglaodblists(species)
pth = fullfile(cdgea, 'assets', 'PanglaoDB', 'celltypes.xlsx');
Tm = readtable(pth, 'FileType', 'spreadsheet', 'Sheet', char(species));
types = string(Tm.CellType);
lists = i_splitlists(string(Tm.PositiveMarkers));
end


function lists = i_splitlists(markerStrings)
lists = cell(1, numel(markerStrings));
for j = 1:numel(markerStrings)
    m = strtrim(split(markerStrings(j), ","));
    lists{j} = m(strlength(m) > 0);
end
end


function [cLtx, types, lists] = i_scorecustom(sce, c, numClusters, Tm)
% Best-scoring list per cluster, by PKG.E_SCORECUSTOMMARKERS (which says
% why it scores against the other clusters rather than on raw counts). A
% cluster that matches no list - best score not above zero - is
% "undetermined".
[score, types, lists] = pkg.e_scorecustommarkers(sce.X, sce.g, c, Tm);
[best, bestType] = max(score, [], 1);
cLtx = repmat("undetermined", numClusters, 1);
isAssigned = best(:) > 0;
cLtx(isAssigned) = types(bestType(isAssigned));
end


function fdr = i_listfdr(sce, c, cLtx, types, lists)
% A hypergeometric check on a marker-list label, as in SAHA's marker-based
% test: do the cluster's own markers overlap the chosen type's marker list
% more than chance would? Both marker-list scorers rank types by a weighted
% expression score, which always has a winner and says nothing about
% whether the winner beats the rest; this is the missing significance. The
% family is every type against every cluster, as in SC_CELLTYPEANNOREF.
queryMarkers = pkg.i_clustermarkers(sce.X, sce.g, c);
P = pkg.e_markeroverlap(queryMarkers, lists, sce.g);
Q = pkg.e_fdr(P);
fdr = nan(numel(cLtx), 1);
for k = 1:numel(cLtx)
    row = types == cLtx(k);
    if any(row)
        fdr(k) = min(Q(row, k));
    end
end
end


function sce = i_recordprovenance(sce, method, detail)
% Without this a saved object cannot say whether its labels came from a
% marker database or from a language model, which are not the same claim.
note = "cell types assigned by " + method + " (" + detail + ") on " + ...
    string(datetime('now', 'Format', 'yyyy-MM-dd HH:mm'));
sce.appendmetainfo(note);
end
