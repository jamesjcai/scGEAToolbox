function [requirerefresh] = callback_SubtypeAnnotation(src, ~)
% CALLBACK_SUBTYPEANNOTATION - Annotate subtypes of an annotated cell type.
%
%   requirerefresh = gui.callback_SubtypeAnnotation(src)
%
% Offers the primary cell types that both the data and the bundled subtype
% marker table know about, then hands each to SC_CSUBTYPEANNO, which isolates
% those cells, re-embeds and re-clusters them on their own, scores the new
% clusters against the subtype markers, and merges the subtype labels back.
%
% Which labels count as a primary type is decided by PKG.I_MATCHPRIMARYTYPE, so
% the menu is available for a gross annotation spelled "T cells_{3}",
% "CD8+ T cells" or "T memory cells", not only for a literal "T cells".
%
% see also: sc_csubtypeanno, pkg.i_matchprimarytype, pkg.i_subtypecandidates,
%           gui.i_updateannotatemenu

requirerefresh = false;

if isa(src, "SingleCellExperiment")
    sce = src;
    FigureHandle = [];
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end

% "undetermined" is what the constructor fills in, not an annotation.
if ~pkg.i_hascelltypelabels(sce)
    gui.myErrordlg(FigureHandle, ['Cells are not annotated yet. Annotate ' ...
        'cell types first, then come back for their subtypes.']);
    return;
end
labels = string(sce.c_cell_type_tx);

% The same question the Annotate menu asks before enabling this item.
[ctypelist, matched, primarytypes] = pkg.i_subtypecandidates(sce);
if isempty(ctypelist)
    gui.myErrordlg(FigureHandle, sprintf(['None of the %d cell type label(s) ' ...
        'in your data is a cell type that cellsubtypes.xlsx has subtype ' ...
        'markers for. Supported: %s.'], numel(unique(labels)), ...
        strjoin(primarytypes, ', ')));
    return;
end

% Say how many cells each choice covers, and under how many different labels:
% the count is what tells the user that "T cells" here also takes in the cells
% annotated as "CD8+ T cells".
items = strings(size(ctypelist));
for k = 1:numel(ctypelist)
    isk = matched == ctypelist(k);
    items(k) = sprintf('%s (%d cells, %d label(s))', ctypelist(k), ...
        sum(isk), numel(unique(labels(isk))));
end

prompt = ['Cells of the selected type are isolated, re-embedded and ' ...
    're-clustered on their own before their subtypes are annotated.'];
if gui.i_isuifig(FigureHandle)
    [indx2, tf2] = gui.myListdlg(FigureHandle, cellstr(items), ...
        'Subtype Annotation', [], true, true, [420, 260], prompt);
else
    [indx2, tf2] = listdlg('PromptString', {prompt}, ...
        'SelectionMode', 'multiple', 'ListString', ...
        cellstr(items), 'ListSize', [420, 260]);
end

if tf2 ~= 1, return; end

celltypetarget_list = ctypelist(indx2);

answer = gui.myQuestdlg(FigureHandle, 'How to label cell type with subtype', ...
    'Choose format', {'Type_{Subtype}','Type (Subtype)','Subtype'}, ...
    'Type_{Subtype}');
if isempty(answer), return; end
switch answer
    case 'Type_{Subtype}'
        formatid = 1;
    case 'Type (Subtype)'
        formatid = 2;
    case 'Subtype'
        formatid = 0;
end

% Re-embedding and re-clustering a subset takes a while, and doing it for
% several types in a row takes several times as long.
fw = gui.myWaitbar(FigureHandle);
for k = 1:length(celltypetarget_list)
    fw = gui.myWaitbar(FigureHandle, fw, false, '', ...
        sprintf('Annotating subtypes of %s...', celltypetarget_list(k)), ...
        (k-1)/length(celltypetarget_list));
    [sce] = sc_csubtypeanno(sce, celltypetarget_list(k), formatid);
end
gui.myWaitbar(FigureHandle, fw);
gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;
end
