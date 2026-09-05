function requirerefresh = callback_AssignCellTypeFromAttrib(src)
% CALLBACK_ASSIGNCELLTYPEFROMATTRIB - Assign cell types from a per-cell column.
%
%   requirerefresh = gui.callback_AssignCellTypeFromAttrib(src)
%
% Picks a per-cell table, then asks which of its columns holds the cell type
% annotation and writes it to sce.c_cell_type_tx. Use this when an imported
% Seurat/RDS or H5AD object carried its annotation under a column name the
% reader did not recognize.
%
% The default source is this dataset's own sce.list_cell_attributes. Unlike the
% View/Export Cell Attribute Table, the standard fields are deliberately left
% out: offering the existing CellType column would let the picker rank it first
% and assign cell type from itself. A table in the base workspace with one row
% per cell, or a delimited text file, can be used instead.
%
% see also: gui.i_pickcelltypecolumn, pkg.i_guesscelltypecol,
%           gui.i_readtablefile, gui.callback_ViewCellAttributeTable, sc_readrdsfile

requirerefresh = false;

[parentfig, sce] = gui.gui_getfigsce(src);
if isempty(sce) || sce.NumCells == 0, return; end

% Restoring a stashed annotation is one of the uses of this callback, so the
% labels it replaces are stashed too and the round trip works in both
% directions. Ask before any of the pickers open.
if ~gui.i_confirmoverwritecelltype(parentfig, sce), return; end

[t, srcname] = in_picksourcetable(sce, parentfig);
if isempty(t), return; end

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

function [t, srcname] = in_picksourcetable(sce, parentfig)
t = [];
srcname = '';
ncells = sce.NumCells;

% Metadata imported from Seurat/RDS lives on the SCE itself as cell attributes,
% so offer those first; the workspace and file options are for tables brought in
% separately.
scetag = '';
tsce = in_attributetable(sce);
if ~isempty(tsce)
    scetag = sprintf('Use this dataset''s own cell attributes (%d columns)', ...
        width(tsce));
end

% The two external sources are named, and their table picked, the same way as
% in GUI.SC_CELLATTRIBEDITOR: one entry each rather than one entry per
% workspace variable. The flat list read better when a qualifying table
% happened to be in the workspace, but it changed shape with whatever was
% there, and when nothing qualified it silently offered nothing at all -
% GUI.I_PICKWORKSPACETABLE says which variables it found and why they did not
% qualify.
loadtag = 'Read a column from a table file (CSV, TSV or TXT)';
wstag = 'Read a column from a table variable in the workspace';
items = {loadtag, wstag};
if ~isempty(scetag)
    items = [{scetag}, items];
end

% The instruction goes in the prompt label, which wraps; the Title stays
% short because the window title bar clips without saying so.
prompt = ['Where are the per-cell columns that hold the cell type ', ...
    'annotation?'];
if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, items, 'Cell Type Source', [], ...
        false, true, [480, 200], prompt);
else
    [indx, tf] = listdlg('PromptString', {prompt}, 'SelectionMode', 'single', ...
        'ListString', items, 'ListSize', [460, 120]);
end
if tf ~= 1 || isempty(indx), return; end

% scetag is '' when the dataset carries no cell attributes, in which case it is
% not in items either, so the case simply never matches.
switch items{indx}
    case scetag
        t = tsce;
        srcname = 'cell attributes';
    case loadtag
        [t, srcname] = gui.i_readtablefile(parentfig, ncells);
    case wstag
        [t, srcname] = gui.i_pickworkspacetable(parentfig, ncells);
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

