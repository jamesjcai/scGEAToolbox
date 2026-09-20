function [requirerefresh] = callback_CollapseCellSubtypes(src, ~)
% CALLBACK_COLLAPSECELLSUBTYPES - Collapse subtype labels back into their parent.
%
%   requirerefresh = gui.callback_CollapseCellSubtypes(src)
%
% Reached from Annotate > Collapse Cell Subtypes back to Cell Type..., next to the two
% items that produce these labels - this is their inverse.
%
% After a subtype run the labels read 'T cells (Regulatory)' or
% 'T cells_{Regulatory}', and after a cluster annotation they read
% 'T cells_{3}'. This offers to drop that last suffix so the cells go back to
% carrying their main cell type.
%
% A label need not carry a suffix to be a subtype. celltypes.xlsx names 37
% subtypes of the primaries cellsubtypes.xlsx covers, so a primary annotation
% run can return 'Plasma cells' or 'Interneurons' - a B cell and a neuron,
% with nothing in either string to say so - and a gross annotation can arrive
% as 'CD8+ T cells'. GUI.I_LABELPARENTS reads those through
% PKG.I_MATCHPRIMARYTYPE, so they are offered here too. Without them this
% item could only undo labels the toolbox had written itself, and a dataset
% annotated straight out of celltypes.xlsx had no way back to its major cell
% types at all.
%
% The proposed merges are shown and confirmed rather than applied outright,
% because only one of the three kinds is a record of something this toolbox
% did. The brace form is written here and nowhere else, so it is preselected.
% The parenthesized form is genuinely ambiguous: 'Endothelial cells (aorta)'
% is a cell type name in celltypes.xlsx, not a subtype of 'Endothelial cells'
% - while 'Endothelial cells (Arterial)', which SC_CSUBTYPEANNO writes for
% that same primary, is one. And a match made on the name alone is a reading
% of what the label means, not a record at all. Both of the latter are listed
% and left unticked.
%
% Only one level comes off per run. 'T cells_{3} (Regulatory)' merges to
% 'T cells_{3}', and running this again takes it to 'T cells'.
%
% A subtype run annotated with the 'Subtype' format keeps no trace of the
% parent in the label. There is nothing to strip there, and the message says
% where the pre-run labels went instead - GUI.CALLBACK_SUBTYPEANNOTATION
% stashes them as a cell attribute.
%
% NOT to be confused with GUI.CALLBACK_MERGECELLSUBTYPES, which is a different
% thing despite the name: that one merges a subtype annotation IN from another
% SingleCellExperiment, in the workspace or in a .mat file. This one only
% rewrites the labels already in this dataset.
%
% see also: gui.i_labelparents, pkg.i_splitsubtypelabel,
%           gui.callback_SubtypeAnnotation, sc_csubtypeanno,
%           gui.callback_AssignCellTypeFromAttrib,
%           gui.callback_MergeCellSubtypes

requirerefresh = false;

if isa(src, "SingleCellExperiment")
    sce = src;
    FigureHandle = [];
else
    [FigureHandle, sce] = gui.gui_getfigsce(src);
end

if ~pkg.i_hascelltypelabels(sce)
    gui.myHelpdlg(FigureHandle, ['Cells are not annotated, so there are no ' ...
        'subtypes to merge.']);
    return;
end

labels = string(sce.c_cell_type_tx);
[ulabels, ~, back] = unique(labels(:));
counts = accumarray(back, 1);
[uparent, uform] = gui.i_labelparents(ulabels);

ismergeable = uparent ~= ulabels;
if ~any(ismergeable)
    gui.myHelpdlg(FigureHandle, ['No cell type label carries a subtype or ' ...
        'subcluster suffix, and none names a subtype of a cell type with ' ...
        'bundled subtype markers, so there is nothing to merge. Labels ' ...
        'written with the "Subtype" format keep no trace of their parent ' ...
        'type; the ' ...
        'labels from before a subtype run are kept as an ''old_cell_type_N'' ' ...
        'cell attribute, and Annotate > Assign Cell Type from Cell Attribute ' ...
        'Table switches back to them.']);
    return;
end

% The name-based rows say so. Without it a user sees 'Plasma cells -> B
% cells' proposed out of a label with no suffix in it and no way to tell why.
idx = find(ismergeable);
note = strings(size(ulabels));
note(uform == "name") = ", by name";
items = arrayfun(@(k) sprintf('%s  ->  %s   (%s%s)', ulabels(k), ...
    uparent(k), pkg.i_plural(counts(k), 'cell'), note(k)), idx, ...
    'UniformOutput', false);

% Preselected: the brace form only. See the note above - a parenthesized name
% may be the cell type's own and a name-based match is a reading rather than a
% record, and merging either is not undone by running this again.
preferred = find(uform(idx) == "brace");

prompt = ['Select the labels to merge into their parent cell type. ' ...
    'Suffixes written as _{...} are ticked: this toolbox writes those. ' ...
    'Suffixes in brackets are left unticked because a name such as ' ...
    '"Endothelial cells (aorta)" is a cell type in its own right, not a ' ...
    'subtype of "Endothelial cells". Rows marked "by name" carry no suffix ' ...
    'at all and are read from the cell type name itself, such as ' ...
    '"Plasma cells" into "B cells"; they are left unticked too.'];
if gui.i_isuifig(FigureHandle)
    [pick, tf] = gui.myListdlg(FigureHandle, items, 'Merge Subtypes', ...
        preferred, true, true, [460, 320], prompt);
else
    [pick, tf] = listdlg('PromptString', {prompt}, ...
        'SelectionMode', 'multiple', 'ListString', items, ...
        'InitialValue', preferred, 'ListSize', [460, 320]);
end
if tf ~= 1 || isempty(pick), return; end

chosen = idx(pick);
newu = ulabels;
newu(chosen) = uparent(chosen);
newlabels = newu(back);
newlabels = reshape(newlabels, size(sce.c_cell_type_tx));

if isequal(string(sce.c_cell_type_tx), newlabels)
    gui.myHelpdlg(FigureHandle, 'No labels changed.');
    return;
end

nbefore = numel(ulabels);
nafter = numel(unique(newlabels));
sce.c_cell_type_tx = newlabels;

gui.myGuidata(FigureHandle, sce, src);
requirerefresh = true;

gui.myHelpdlg(FigureHandle, sprintf( ...
    '%s merged. The cells now carry %s, down from %d.', ...
    pkg.i_plural(numel(chosen), 'label'), ...
    pkg.i_plural(nafter, 'cell type'), nbefore));
end
