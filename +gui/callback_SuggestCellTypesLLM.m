function callback_SuggestCellTypesLLM(src, ~)
%CALLBACK_SUGGESTCELLTYPESLLM Show a language model's cell type guesses, save nothing.
%
%   gui.callback_SuggestCellTypesLLM(app)
%
% Annotate > Suggest Cell Types with LLM. Runs LLM.E_CELLTYPEANNO on the
% clusters and shows its answers next to the current labels, for the user to
% weigh. Nothing is written to the SCE: no labels, no stash, no undo step.
%
% Suggest-only by design. LLM annotation as a labelling method was removed
% from this menu (commit 9491d20) because a model's answer is not a sound
% basis for labelling cells: it can name a cluster wrongly with full vote
% agreement, as on example_data/new_example_sce.mat, where a cluster with
% myeloid markers came back "T cells" three times out of three. As a second
% opinion beside the marker-based labels it is still useful. To apply its
% answers, use sc_annotatecells(sce, Method="llm") on the command line.
%
% see also: llm.e_celltypeanno, sc_annotatecells, gui.i_setllmmodel

[FigureHandle, sce] = gui.gui_getfigsce(src);
if isempty(sce) || sce.NumCells < 1, return; end

% The method labels clusters. SingleCellExperiment seeds c_cluster_id to
% all-ones, so a single group means not clustered.
if isempty(sce.c_cluster_id) || isscalar(unique(string(sce.c_cluster_id)))
    gui.myWarndlg(FigureHandle, ['The language model is asked about ' ...
        'each cluster''s marker genes, and the cells have not been ' ...
        'clustered. Cluster them first (Cluster > Cluster Cells).']);
    return;
end

% Configured provider. Offer to set one rather than failing on it.
try
    [provider, model] = llm.i_llmsettings();
catch
    answer = gui.myQuestdlg(FigureHandle, ['No LLM provider is set. ' ...
        'Select one now?']);
    if ~strcmp(answer, 'Yes'), return; end
    if ~gui.i_setllmmodel(src), return; end
    [provider, model] = llm.i_llmsettings();
end

speciestag = gui.i_selectspecies(2, false, FigureHandle);
if isempty(speciestag), return; end

% A named tissue changes the answer a lot. For the catalogued ones the model
% is also given that tissue's cell-type list to choose from.
tissues = llm.i_celltissues();
choices = [tissues(:); "Other (type it in)..."; "Not specified"];
[indx, tf] = gui.myListdlg(FigureHandle, choices, ...
    'Tissue or organ (a listed one adds a reference cell-type list):', ...
    [], false);
if tf ~= 1 || isempty(indx), return; end
switch choices(indx(1))
    case "Not specified"
        tissue = "";
    case "Other (type it in)..."
        answer = gui.myInputdlg({'Tissue or organ:'}, ...
            'Suggest Cell Types with LLM', {''}, FigureHandle);
        if isempty(answer), return; end
        tissue = strtrim(string(answer{1}));
    otherwise
        tissue = choices(indx(1));
end

fw = gui.myWaitbar(FigureHandle, [], false, sprintf(['Asking %s / %s ' ...
    'about each cluster (several minutes)...'], provider, model));
try
    [~, ~, T] = llm.e_celltypeanno(sce, Species=string(speciestag), ...
        Tissue=tissue, Provider=provider, Model=model, Verbose=false);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

T = in_addcontext(T, sce);
gui.myHelpdlg(FigureHandle, sprintf(['Suggestions for %d clusters from ' ...
    '%s / %s. Nothing was saved.\n\nAgreement is how often the model ' ...
    'gave the same answer when asked again, not how often it is right. ' ...
    'Check each suggestion against the Markers column.'], ...
    height(T), provider, model));
gui.TableViewerApp(T, FigureHandle, "LLMCellTypeSuggestions");
end


function T = in_addcontext(T, sce)
% Cluster size, and the cluster's current majority label when the cells are
% annotated, placed after CellType so the two read side by side.
cid = string(sce.c_cluster_id(:));
cluster = string(T.Cluster);
numCells = zeros(height(T), 1);
current = strings(height(T), 1);
hasLabels = pkg.i_hascelltypelabels(sce);
for k = 1:height(T)
    inCluster = cid == cluster(k);
    numCells(k) = sum(inCluster);
    if hasLabels && numCells(k) > 0
        labels = string(sce.c_cell_type_tx(inCluster));
        [u, ~, j] = unique(labels);
        [~, top] = max(accumarray(j, 1));
        current(k) = u(top);
    end
end
T = addvars(T, numCells, 'After', 'Cluster', 'NewVariableNames', 'NumCells');
if hasLabels
    T = addvars(T, current, 'After', 'CellType', ...
        'NewVariableNames', 'CurrentLabel');
end
end
