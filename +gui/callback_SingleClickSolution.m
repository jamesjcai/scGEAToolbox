function [requirerefresh, speciestag] = callback_SingleClickSolution(src, ~)

requirerefresh = false;
speciestag = [];

[FigureHandle, sce] = gui.gui_getfigsce(src);

if isa(src, 'matlab.apps.AppBase')
    speciestag = src.speciestag;
end

if ~isprop(sce, 'c_cell_type_tx')
    disp('The sce object does not have the property c_cell_type_tx.');
    return;
end
if ~all(sce.c_cell_type_tx == "undetermined")
    if ~strcmp(gui.myQuestdlg(FigureHandle, ...
            "Your data has been embedded and annotated. " + ...
            "Single Click Solution will re-embed and " + ...
            "annotate cells. Current embedding and " + ...
            "annotation will be overwritten. Continue?"), 'Yes')
        return;
    end
else
    if ~gui.gui_showrefinfo('Single Click Solution', FigureHandle), return; end
end

hasDuplicates = numel(unique(sce.g)) < numel(sce.g);
if hasDuplicates
    [~, sce] = gui.gui_rmdugenes(sce, FigureHandle);
end

speciestag = gui.i_selectspecies(2, false, FigureHandle, speciestag);

if isempty(speciestag) || strlength(speciestag) == 0, return; end

% The gene-set item sits with the embeddings it modifies rather than at the
% end of the list. It is a checkbox and not GUI.I_GETHVGNUM's full four-way
% choice because this is a preset pipeline - another modal dialog would
% undercut "single click" - and it is off by default, which is the
% behaviour this callback has always had.
%
% The items are indexed by name. This block used to walk PROMPT with a
% COUNT that was incremented between the steps, which silently ties every
% later step's identity to the length of the list above it: inserting the
% gene-set item shifted cell cycle and potency onto each other's labels.
item.tsne = 1;
item.umap = 2;
item.phate = 3;
item.markers = 4;
item.cellcycle = 5;
item.potency = 6;

prompt = {
    'tSNE Embedding?', ...
    'Add UMAP Embedding?', ...
    'Add PHATE Embedding?', ...
    'Include PanglaoDB Markers in Embedding?', ...
    'Estimate Cell Cycles?', ...
    'Estimate Differentiation Potency of Cells?'};
assert(numel(prompt) == numel(fieldnames(item)), ...
    'CALLBACK_SINGLECLICKSOLUTION: every PROMPT item needs an ITEM index.');

answer = gui.myChecklistdlg(FigureHandle, prompt, ...
    'Title', 'Select Items', 'DefaultSelection', item.tsne);

if isempty(answer)
    return;
end

if ~ismember(prompt{item.tsne}, answer)
    gui.myErrordlg(FigureHandle, 'tSNE Embedding has to be included.', '');
    return;
end

% The markers go into the clustering's principal components as well as
% into the embeddings, or ticking this would change only the display.
markergenes = [];
if ismember(prompt{item.markers}, answer)
    genemode = "hvg+markers";
    markergenes = pkg.i_getmarkerwhitelist(sce.g, sce.X);
else
    genemode = "hvg";
end

fw = gui.myWaitbar(FigureHandle);
gui.myWaitbar(FigureHandle, fw, false, '', 'Basic QC Filtering...', 1/8);
sce = sce.qcfilter;

gui.myWaitbar(FigureHandle, fw, false, '', 'Embedding cells using tSNE...', 2/8);
try
    sce = sce.embedcells('tsne3d', true, genemode, 3);
catch ME
    gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

if ismember(prompt{item.umap}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', 'Embedding cells using UMAP...', 2/8);
    try
        sce = sce.embedcells('umap3d', true, genemode, 3);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
end

if ismember(prompt{item.phate}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', 'Embedding cells using PHATE...', 2/8);
    try
        sce = sce.embedcells('phate3d', true, genemode, 3);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
end

% Louvain on an SNN graph of the principal components, as Seurat's
% FindNeighbors/FindClusters do, with the resolution (0.8) rather than a
% preset k deciding the number of clusters. This used to be k-means on the
% 3-D t-SNE, which inherits t-SNE's distortions and its random start. On
% the mouse pancreas example (8,260 cells, 11 labelled types) the switch
% took the ARI against the labels from 0.11 (k = 80) to 0.53 (20
% clusters), the ARI between reruns from 0.69 to 1.00, and annotation
% accuracy from 0.89 to 0.90. The t-SNE is now for display only.
% It is stored as LOUVAINPC, not LOUVAIN, which is Louvain on the embedding.
gui.myWaitbar(FigureHandle, fw, false, '', ...
    'Clustering cells using Louvain on principal components...', 3/8);
sce = sce.clustercells([], 'louvainpc', true, [], Genes=markergenes);
gui.myWaitbar(FigureHandle, fw, false, '', 'Annotating cell types using PanglaoDB...', 4/8);
sce = sce.assigncelltype(speciestag, false);

if ismember(prompt{item.cellcycle}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', 'Estimate cell cycles...', 5/8);
    sce = sce.estimatecellcycle;
end

if ismember(prompt{item.potency}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', ...
        'Estimate differentiation potency of cells...', 6/8);
    sce = sce.estimatepotency(speciestag);
end

gui.myWaitbar(FigureHandle, fw, false, '', '', 7/8);
gui.myGuidata(FigureHandle, sce, src);
gui.myWaitbar(FigureHandle, fw);

requirerefresh = true;
end
