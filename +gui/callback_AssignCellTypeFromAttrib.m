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

[t, srcname] = in_picksourcetable(sce, parentfig);
if isempty(t), return; end

% Score the columns before anything is put on screen. When nothing qualifies the
% warning must be the only window that opens: myWarndlg draws inside parentfig
% rather than as its own window, so a table viewer opened first would hide it.
[~, ~, ranked] = pkg.i_guesscelltypecol(t, sce.NumCells);
if isempty(ranked) || height(ranked) == 0
    gui.myWarndlg(parentfig, in_noqualifyreason(t, srcname, sce.NumCells));
    return;
end

% Show the candidate columns alongside the picker so their values can be
% inspected while choosing. Only the ranked columns are shown, so the table
% matches the choices on offer instead of also listing columns the picker has
% already ruled out. drawnow realizes the window now, otherwise its deferred
% paint flushes while the picker waits and raises it over the picker.
viewerfig = gui.TableViewerApp(t(:, ranked.Name), parentfig, "CellTypeCandidates");
drawnow;

[ctype, colname] = gui.i_pickcelltypecolumn(t, sce, parentfig, [], ranked);

% Close the viewer before any message dialog, which is drawn inside parentfig
% and would otherwise sit behind the viewer window.
if pkg.i_isvalid(viewerfig), delete(viewerfig); end
if isempty(ctype), return; end

sce.c_cell_type_tx = ctype;
gui.myGuidata(parentfig, sce, src);
requirerefresh = true;

gui.myHelpdlg(parentfig, sprintf( ...
    'Cell type assigned from %s, column "%s" (%d types).', ...
    srcname, colname, numel(unique(ctype))));
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
    scetag = sprintf('Cell attributes on this dataset (%d columns)', width(tsce));
end

names = {};
if ~(ismcc || isdeployed)
    a = evalin('base', 'whos');
    if ~isempty(a)
        istbl = strcmp({a.class}, 'table');
        nrow = cellfun(@(x) x(1), {a.size});
        names = {a(istbl & nrow == ncells).name};
    end
end

loadtag = 'Load Table from File...';
items = [names, {loadtag}];
if ~isempty(scetag)
    items = [{scetag}, items];
end

prompt = 'Select a source of per-cell columns:';
if gui.i_isuifig(parentfig)
    [indx, tf] = gui.myListdlg(parentfig, items, prompt, [], false);
else
    [indx, tf] = listdlg('PromptString', {prompt}, ...
        'SelectionMode', 'single', 'ListString', items, 'ListSize', [300, 300]);
end
if tf ~= 1 || isempty(indx), return; end

if ~isempty(scetag) && strcmp(items{indx}, scetag)
    t = tsce;
    srcname = 'cell attributes';
elseif strcmp(items{indx}, loadtag)
    [t, srcname] = gui.i_readtablefile(parentfig, ncells);
    if isempty(t), return; end
else
    t = evalin('base', items{indx});
    srcname = items{indx};
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
