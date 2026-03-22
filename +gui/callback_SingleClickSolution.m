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

prompt = {
    'tSNE Embedding?', ...
    'Add UMAP Embedding?', ...
    'Add PHATE Embedding?', ...
    'Estimate Cell Cycles?', ...
    'Estimate Differentiation Potency of Cells?'};

answer = gui.myChecklistdlg(FigureHandle, prompt, ...
    'Title', 'Select Items', 'DefaultSelection', 1);

if isempty(answer)
    return;
end

count = 1;
if ~ismember(prompt{count}, answer)
    gui.myErrordlg(FigureHandle, 'tSNE Embedding has to be included.', '');
    return;
end

fw = gui.myWaitbar(FigureHandle);
gui.myWaitbar(FigureHandle, fw, false, '', 'Basic QC Filtering...', 1/8);
sce = sce.qcfilter;

gui.myWaitbar(FigureHandle, fw, false, '', 'Embedding cells using tSNE...', 2/8);
try
    sce = sce.embedcells('tsne3d', true, true, 3);
catch ME
    gui.myWaitbar(FigureHandle, fw);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end

count = count + 1;
if ismember(prompt{count}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', 'Embedding cells using UMAP...', 2/8);
    try
        sce = sce.embedcells('umap3d', true, true, 3);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
end

count = count + 1;
if ismember(prompt{count}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', 'Embedding cells using PHATE...', 2/8);
    try
        sce = sce.embedcells('phate3d', true, true, 3);
    catch ME
        gui.myWaitbar(FigureHandle, fw);
        gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
        return;
    end
end

gui.myWaitbar(FigureHandle, fw, false, '', 'Clustering cells using K-means...', 3/8);
sce = sce.clustercells([], [], true);
gui.myWaitbar(FigureHandle, fw, false, '', 'Annotating cell types using PanglaoDB...', 4/8);
sce = sce.assigncelltype(speciestag, false);

count = count + 1;
if ismember(prompt{count}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', 'Estimate cell cycles...', 5/8);
    sce = sce.estimatecellcycle;
end

count = count + 1;
if ismember(prompt{count}, answer)
    gui.myWaitbar(FigureHandle, fw, false, '', ...
        'Estimate differentiation potency of cells...', 6/8);
    sce = sce.estimatepotency(speciestag);
end

gui.myWaitbar(FigureHandle, fw, false, '', '', 7/8);
gui.myGuidata(FigureHandle, sce, src);
gui.myWaitbar(FigureHandle, fw);

requirerefresh = true;
end
