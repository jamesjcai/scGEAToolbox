function requirerefresh = callback_AssignCellTypeFromAttrib(src)
% CALLBACK_ASSIGNCELLTYPEFROMATTRIB - Assign cell types from a per-cell column.
%
%   requirerefresh = gui.callback_AssignCellTypeFromAttrib(src)
%
% Asks which column of this dataset's own cell attributes
% (sce.list_cell_attributes) holds the cell type annotation and writes it to
% sce.c_cell_type_tx. Use this when an imported Seurat/RDS or H5AD object
% carried its annotation under a column name the reader did not recognize.
%
% Unlike the View/Export Cell Attribute Table, the standard fields are
% deliberately left out: offering the existing CellType column would let the
% picker rank it first and assign cell type from itself. A column held in a
% table file or workspace variable is brought in first with Edit > Add/Edit
% Cell Attributes, which reads both.
%
% see also: gui.i_pickcelltypecolumn, pkg.i_guesscelltypecol,
%           gui.sc_cellattribeditor, gui.callback_ViewCellAttributeTable, sc_readrdsfile

requirerefresh = false;

[parentfig, sce] = gui.gui_getfigsce(src);
if isempty(sce) || sce.NumCells == 0, return; end

% The menu item names its source, Annotate > Import Annotations > from a Cell
% Attribute Column, so there is no source picker: an empty table is reported
% here, before anything else is asked, rather than falling through to a file
% or workspace prompt the user did not ask for.
t = in_attributetable(sce);
srcname = 'the cell attribute table';
if width(t) == 0
    gui.myWarndlg(parentfig, in_noattributesreason(sce));
    return;
end

% Restoring a stashed annotation is one of the uses of this callback, so the
% labels it replaces are stashed too and the round trip works in both
% directions. Ask before the column picker opens.
if ~gui.i_confirmoverwritecelltype(parentfig, sce), return; end

% Score the columns before anything is put on screen, so that when nothing
% qualifies the warning is the only window that opens.
[~, ~, ranked] = pkg.i_guesscelltypecol(t, sce.NumCells);
if isempty(ranked) || height(ranked) == 0
    gui.myWarndlg(parentfig, in_noqualifyreason(t, srcname, sce.NumCells));
    return;
end

% The picker alone, no table viewer alongside it. A viewer of the candidate
% columns was shown here so their values could be inspected while choosing,
% but it is a separate window while the picker is drawn inside parentfig, so
% the two could not be placed or stacked without one obscuring the other. The
% picker already lists each candidate's distinct-value count and example
% values, which is what the viewer was there to supply.
[ctype, colname] = gui.i_pickcelltypecolumn(t, sce, parentfig, [], ranked);
if isempty(ctype), return; end

stashname = pkg.i_stashcelltypehistory(sce);
sce.c_cell_type_tx = ctype;
gui.myGuidata(parentfig, sce, src);
requirerefresh = true;

% Recolor the main plot by the labels just assigned. The menu handler in
% scgeatoolApp.mlapp ignores the returned REQUIREREFRESH, so without this the
% points keep whatever grouping they carried before and the assignment looks
% as though it did nothing. GUI.CALLBACK_ANNOTATECELLS does the same.
if isa(src, 'matlab.apps.AppBase')
    [src.c, src.cL] = findgroups(string(src.sce.c_cell_type_tx));
    src.sce.c = src.c;
    src.in_RefreshAll(true, false);
    src.ix_labelclusters(true);
end

msg = sprintf('Cell type assigned from %s, column "%s" (%d types).', ...
    srcname, colname, numel(unique(ctype)));
gui.myHelpdlg(parentfig, msg + gui.i_stashnotice(stashname));
end

function msg = in_noqualifyreason(t, srcname, ncells)
% Say why nothing qualified. The columns the picker can score are a subset of
% what the table holds, so "no column looks like a cell type" on its own reads
% as a bug when the table plainly has columns in it.

if width(t) == 0
    msg = sprintf('There are no columns in %s to choose from.', srcname);
    return;
end

istext = false(1, width(t));
for k = 1:width(t)
    v = t.(k);
    istext(k) = (isstring(v) || iscellstr(v) || iscategorical(v)) && size(v, 2) == 1;
end

if ~any(istext)
    msg = sprintf(['None of the %d columns of %s holds text labels. A cell ' ...
        'type column must contain text, not numbers.'], width(t), srcname);
else
    % Mirrors the cap applied in pkg.i_guesscelltypecol.
    maxtypes = max(2, min(200, ceil(ncells/5)));
    msg = sprintf(['None of the %d text columns of %s looks like a cell type ' ...
        'annotation. A column qualifies when it has between 2 and %d distinct ' ...
        'labels and those labels are not all numbers.'], ...
        sum(istext), srcname, maxtypes);
end
end

function msg = in_noattributesreason(sce)
% Say why there is nothing to pick from, and where a column can come from.

hint = ['Add one with Edit > Add/Edit Cell Attributes, which can read it ' ...
    'from a table file or a workspace variable, then try again.'];
if isempty(sce.list_cell_attributes)
    msg = ['This dataset has no cell attributes beyond the standard fields. ' hint];
elseif numel(sce.list_cell_attributes) == 2
    msg = sprintf(['The only cell attribute of this dataset does not have ' ...
        'one value per cell (%d cells). %s'], sce.NumCells, hint);
else
    msg = sprintf(['None of the %d cell attributes of this dataset has one ' ...
        'value per cell (%d cells). %s'], ...
        numel(sce.list_cell_attributes)/2, sce.NumCells, hint);
end
end

function t = in_attributetable(sce)
% Build a table from list_cell_attributes only. The standard SCE fields are
% deliberately left out: including the existing CellType column would let the
% picker score it top and assign cell type from itself.

t = table();
if isempty(sce.list_cell_attributes), return; end

names = string(sce.list_cell_attributes(1:2:end));
vals = sce.list_cell_attributes(2:2:end);
isok = cellfun(@(v) numel(v) == sce.NumCells, vals);
if ~any(isok), return; end

names = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names(isok)));
vals = cellfun(@(v) v(:), vals(isok), 'UniformOutput', false);
t = table(vals{:}, 'VariableNames', names);
end

