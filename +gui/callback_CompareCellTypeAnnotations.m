function callback_CompareCellTypeAnnotations(src, ~)
%CALLBACK_COMPARECELLTYPEANNOTATIONS Show the cell type annotations side by side.
%
%   Every annotation channel stashes the labels it replaces as a numbered
%   'old_cell_type_N' cell attribute (PKG.I_STASHCELLTYPEHISTORY), so a dataset
%   annotated more than once carries a history. Keeping it is only half the
%   job: this is how it gets looked at.
%
%   Three views, because the useful question differs. The plots answer "where
%   do these disagree" spatially; the table answers "what did each method call
%   THIS cell"; the cross-tabulation answers "which types got split or merged",
%   which neither of the other two shows.
%
%   See also PKG.I_CELLTYPEHISTORY, PKG.I_STASHCELLTYPEHISTORY,
%   GUI.I_ADDCOMPAREANNOTMENU.

[parentfig, sce] = gui.gui_getfigsce(src);
if isempty(sce) || sce.NumCells == 0, return; end

[names, labels] = pkg.i_celltypehistory(sce);
if numel(names) < 2
    gui.myHelpdlg(parentfig, ['This dataset carries only one cell type ', ...
        'annotation, so there is nothing to compare. Annotating again keeps ', ...
        'the current labels as an ''old_cell_type_N'' cell attribute, and ', ...
        'they then show up here.']);
    return;
end

viewitems = { ...
    'Side-by-side plots on the cell embedding', ...
    'Table with one row per cell', ...
    'Cross-tabulate two annotations (what got split or merged)'};
prompt = sprintf(['This dataset carries %d cell type annotations. ', ...
    'How should they be shown?'], numel(names));
[indx, tf] = gui.myListdlg(parentfig, viewitems, 'Cell Type Annotations', ...
    viewitems{1}, false, true, [480, 180], prompt);
if tf ~= 1 || isempty(indx), return; end

switch indx
    case 1
        in_plotpanels(sce, names, labels, parentfig);
    case 2
        in_showtable(sce, names, labels, parentfig);
    case 3
        in_crosstab(names, labels, parentfig);
end
end


function in_plotpanels(sce, names, labels, parentfig)
% One embedding panel per annotation, on a shared set of axes limits.
%
% A shared limit is what makes the panels comparable: with each panel free to
% autoscale, the same cell sits at a different place in each one and the eye
% cannot follow it.
%
% GSCATTER rather than GUI.I_GSCATTER3, which maps the labels through
% findgroups and so cannot label a legend with the type names. An earlier
% version wrote the names at each group's centroid instead, which reads well
% for well-separated groups and turns to overlapping mush as soon as two types
% share a blob - exactly the case this view exists to show. A legend never
% collides. Past in_maxlegendtypes() types a legend is taller than the panel,
% so it is dropped and the type count in the title carries what is left.

if size(sce.s, 1) ~= sce.NumCells || size(sce.s, 2) < 2
    gui.myWarndlg(parentfig, ['The cells have no 2-D embedding to plot on. ', ...
        'Run an embedding first, or use the table view.']);
    return;
end

n = numel(names);
ncol = ceil(sqrt(n));
nrow = ceil(n/ncol);
f = figure('Name', 'Cell Type Annotations', 'NumberTitle', 'off');
tl = tiledlayout(f, nrow, ncol, 'TileSpacing', 'compact', 'Padding', 'compact');

s = sce.s;
xl = [min(s(:,1)), max(s(:,1))];
yl = [min(s(:,2)), max(s(:,2))];
pad = 0.03*[-1 1];
xl = xl + pad*diff(xl);
yl = yl + pad*diff(yl);

for k = 1:n
    ax = nexttile(tl);
    lbl = labels{k};
    ntypes = numel(unique(lbl));
    gscatter(ax, s(:,1), s(:,2), lbl, [], '.', 8);
    if ntypes <= in_maxlegendtypes()
        lg = legend(ax, 'Location', 'eastoutside');
        lg.Interpreter = 'none';
        lg.FontSize = 7;
        lg.Box = 'off';
    else
        legend(ax, 'off');
    end
    title(ax, sprintf('%s  (%d %s)', names(k), ntypes, ...
        in_plural('type', ntypes)), 'Interpreter', 'none');
    xlim(ax, xl); ylim(ax, yl);
    ax.XTick = []; ax.YTick = [];
    xlabel(ax, ''); ylabel(ax, '');
end
end


function n = in_maxlegendtypes()
n = 14;
end


function s = in_plural(word, n)
if n == 1
    s = word;
else
    s = [word 's'];
end
end

function in_showtable(sce, names, labels, parentfig)
% One row per cell, one column per annotation, cell barcode first.

vars = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(names));
t = table(string(sce.c_cell_id(:)), 'VariableNames', {'CellID'});
for k = 1:numel(names)
    t.(vars(k)) = labels{k};
end

% A column saying whether the annotations agree on a cell is the reason to
% look at this table at all, so it is filled in rather than left to the eye.
same = true(sce.NumCells, 1);
for k = 2:numel(labels)
    same = same & (labels{k} == labels{1});
end
t.Agree = same;

gui.TableViewerApp(t, parentfig, 'CellTypeAnnotations');
end


function in_crosstab(names, labels, parentfig)
% Cross-tabulate two annotations: rows one, columns the other, counts inside.

[indx, tf] = gui.myListdlg(parentfig, cellstr(names), 'Cross-tabulate', ...
    [], true, true, [420, 260], ['Select exactly two annotations. Rows will ', ...
    'be the first, columns the second.']);
if tf ~= 1 || numel(indx) ~= 2
    if tf == 1
        gui.myWarndlg(parentfig, sprintf(['Select exactly two annotations ', ...
            'to cross-tabulate; %d were selected.'], numel(indx)));
    end
    return;
end

a = labels{indx(1)};
b = labels{indx(2)};
[ga, la] = findgroups(a);
[gb, lb] = findgroups(b);
M = accumarray([ga(:), gb(:)], 1, [numel(la), numel(lb)]);

colnames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(lb));
t = array2table(M, 'VariableNames', colnames);
t = addvars(t, la(:), 'Before', 1, 'NewVariableNames', {'RowLabel'});
gui.TableViewerApp(t, parentfig, 'CellTypeCrosstab');

% Exact-label agreement is only meaningful when the two annotations share a
% vocabulary, which two different methods often do not, so say which it is
% rather than reporting a number that looks worse than the result is.
shared = intersect(la, lb);
if isempty(shared)
    msg = sprintf(['"%s" and "%s" share no label names, so only the ', ...
        'cross-tabulation is meaningful - a per-cell agreement rate would ', ...
        'read as 0%% however well the groupings line up.'], ...
        names(indx(1)), names(indx(2)));
else
    msg = sprintf(['"%s" and "%s" give the same label to %.1f%% of cells ', ...
        '(%d of %d), over %d shared label name(s).'], ...
        names(indx(1)), names(indx(2)), 100*mean(a == b), sum(a == b), ...
        numel(a), numel(shared));
end
gui.myHelpdlg(parentfig, msg);
end
