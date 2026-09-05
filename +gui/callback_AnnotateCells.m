function [needupdatesce] = callback_AnnotateCells(src, ~)
%CALLBACK_ANNOTATECELLS Assign cell types, choosing the method in one place.
%
%   One menu item covering every annotation method the toolbox has, instead
%   of one item per method where the choice is made by which menu you open.
%   Marker-gene-based annotation has its own, "Annotate Cell Types Using
%   PanglaoDB Marker Genes" and "...Customized Marker Genes...", which still
%   work; this does not replace them. SCimilarity and panhumanpy delegate to
%   GUI.CALLBACK_RUNSCIMILARITY and GUI.CALLBACK_RUNPANHUMANPY - the same
%   code External > Atlas-Based Cell Type Annotation uses - rather than
%   duplicating a thinner path through SC_ANNOTATECELLS.
%
%   See also SC_ANNOTATECELLS, GUI.CALLBACK_RUNSCIMILARITY,
%   GUI.CALLBACK_RUNPANHUMANPY, GUI.I_ADDANNOTATEMENU.

needupdatesce = false;
[FigureHandle, sce] = gui.gui_getfigsce(src);
if isempty(sce) || sce.NumCells < 1, return; end

% ---- method -------------------------------------------------------------
methodtag = ["llm", "consensus", "scimilarity", "panhumanpy"];
descr = [ ...
    "Language model reading each cluster's marker genes"
    "Consensus of markers + language model, disagreements flagged"
    "SCimilarity reference model (needs a downloaded model)"
    "panhumanpy reference model (needs Python)"];

[indx, tf] = gui.myListdlg(FigureHandle, descr, ...
    'Annotate cell types using:', [], false);
if tf ~= 1 || isempty(indx), return; end
method = methodtag(indx(1));

% ---- reference-model methods: reuse their full-featured callbacks --------
% These work per cell (no clustering needed) and already handle model/
% target-cell-type selection, the memory check, the working folder, and
% (like every path here) confirming before they overwrite existing labels.
if ismember(method, ["scimilarity", "panhumanpy"])
    switch method
        case "scimilarity"
            needupdatesce = gui.callback_RunSCimilarity(src, []);
        case "panhumanpy"
            needupdatesce = gui.callback_RunPanhumanpy(src, []);
    end
    if needupdatesce
        app = src;
        [app.c, app.cL] = findgroups(string(app.sce.c_cell_type_tx));
        app.sce.c = app.c;
        in_RefreshAll(app, true, false);
        ix_labelclusters(app, true);
    end
    return;
end

% ---- warn before overwriting existing labels (llm, consensus) -----------
if ~gui.i_confirmoverwritecelltype(FigureHandle, sce), return; end

% ---- what that method needs --------------------------------------------
opts = struct('Method', method, 'Verbose', true);

% These label clusters. Say so now rather than after the user has answered
% two more dialogs and waited for a model to run.
if isempty(sce.c_cluster_id) || isscalar(unique(sce.c_cluster_id))
    gui.myWarndlg(FigureHandle, ['This method labels clusters, and ', ...
        'the data has none. Cluster the cells first, or pick a ', ...
        'reference-model method, which works per cell.']);
    return;
end
speciestag = gui.i_selectspecies(2, false, FigureHandle);
if isempty(speciestag), return; end
opts.Species = string(speciestag);

% Naming the tissue changes the answer a lot: the same marker set is far
% less ambiguous once the tissue is known. For the catalogued tissues it
% does more than that - the model is given that tissue's full cell-type
% list to choose from, instead of naming one from memory.
tissues = llm.i_celltissues();
choices = [tissues(:); "Other (type it in)..."; "Not specified"];
[indx, tf] = gui.myListdlg(FigureHandle, choices, ...
    'Tissue or organ (a listed one adds a reference cell-type list):', ...
    [], false);
if tf ~= 1 || isempty(indx), return; end

switch choices(indx(1))
    case "Not specified"
        opts.Tissue = "";
    case "Other (type it in)..."
        answer = gui.myInputdlg({'Tissue or organ:'}, ...
            'Annotate Cell Types', {''}, FigureHandle);
        if isempty(answer), return; end
        opts.Tissue = string(answer{1});
    otherwise
        opts.Tissue = choices(indx(1));
end

% ---- run ----------------------------------------------------------------
fw = gui.myWaitbar(FigureHandle, [], false, ...
    sprintf('Annotating cell types using %s...', method));
try
    args = namedargs2cell(opts);
    [sce, T, stashname] = sc_annotatecells(sce, args{:});
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

if isempty(T) || height(T) == 0
    gui.myHelpdlg(FigureHandle, 'No cell types were assigned.');
    return;
end

gui.myGuidata(FigureHandle, sce, src);
needupdatesce = true;

% ---- report -------------------------------------------------------------
% Uncertainty is the most useful thing consensus produces, so put it in the
% dialog rather than leaving it in a table the user may never open.
msg = sprintf('%d cell type(s) assigned to %d cells using "%s".', ...
    numel(unique(T.CellType)), sce.NumCells, method);
if any(~isnan(T.Confidence))
    nunsure = sum(T.Confidence < 1);
    if nunsure == 0
        msg = sprintf('%s\nEvery cluster was called with full agreement.', msg);
    else
        msg = sprintf(['%s\n%d of %d cluster(s) were uncertain or had ', ...
            'methods disagreeing - see the Confidence column.'], ...
            msg, nunsure, height(T));
    end
end
gui.myHelpdlg(FigureHandle, msg + gui.i_stashnotice(stashname));

gui.i_exporttable(T, true, 'Tcelltypeanno', 'CellTypeAnnotation', ...
    [], [], FigureHandle);

app = src;
[app.c, app.cL] = findgroups(string(app.sce.c_cell_type_tx));
app.sce.c = app.c;
in_RefreshAll(app, true, false);
ix_labelclusters(app, true);
end
