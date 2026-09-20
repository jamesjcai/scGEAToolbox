function callback_CompareCellTypeAnnotations(src, ~)
%CALLBACK_COMPARECELLTYPEANNOTATIONS Show the cell type annotations side by side.
%
%   Every annotation channel stashes the labels it replaces as a numbered
%   'old_cell_type_N' cell attribute (PKG.I_STASHCELLTYPEHISTORY), so a dataset
%   annotated more than once carries a history. Keeping it is only half the
%   job: this is how it gets looked at.
%
%   Four views, because the useful question differs. The plots answer "where
%   do these disagree" spatially; the table answers "what did each method call
%   THIS cell"; the cross-tabulation and the flow diagram both answer "which
%   types got split or merged", the first exactly and the second at a glance.
%
%   See also PKG.I_CELLTYPEHISTORY, PKG.I_STASHCELLTYPEHISTORY,
%   GUI.I_UPDATEANNOTATEMENU.

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
    'Side-by-side plots on the cell embedding (linked by brushing)', ...
    'Table with one row per cell', ...
    'Cross-tabulate two annotations (what got split or merged)', ...
    'Sankey flow diagram between two annotations'};
prompt = sprintf(['This dataset carries %d cell type annotations. ', ...
    'How should they be shown?'], numel(names));
[indx, tf] = gui.myListdlg(parentfig, viewitems, 'Cell Type Annotations', ...
    viewitems{1}, false, true, [480, 200], prompt);
if tf ~= 1 || isempty(indx), return; end

switch indx
    case 1
        in_plotpanels(sce, names, labels, parentfig);
    case 2
        in_showtable(sce, names, labels, parentfig);
    case 3
        in_crosstab(names, labels, parentfig);
    case 4
        in_sankey(names, labels, parentfig);
end
end


function pick = in_picktwo(names, parentfig, prompt)
% The two views that compare a PAIR share this picker, so they cannot drift
% into disagreeing about what "exactly two" means or in which order.
pick = [];
[indx, tf] = gui.myListdlg(parentfig, cellstr(names), 'Select two', ...
    [], true, true, [420, 260], prompt);
if tf ~= 1 || numel(indx) ~= 2
    if tf == 1
        gui.myWarndlg(parentfig, sprintf(['Select exactly two ', ...
            'annotations; %d were selected.'], numel(indx)));
    end
    return;
end
pick = indx;
end


function in_sankey(names, labels, parentfig)
% Ribbons from each type in the first annotation to each type in the second.
pick = in_picktwo(names, parentfig, ['Select exactly two annotations. The ', ...
    'flow runs from the first to the second.']);
if isempty(pick), return; end

gui.i_alluvialview(labels{pick(1)}, labels{pick(2)}, names(pick(1)), ...
    names(pick(2)), parentfig, 'Cell Type Annotation Flow');
end


function in_plotpanels(sce, names, labels, parentfig)
% One embedding panel per annotation, drawn by the Multi-Grouping View.
%
% GUI.I_MULTIGROUPVIEW is the same figure that Group > Multi-Grouping View
% opens, handed this dataset's annotation history instead of a picker's
% selection. Writing a second panel drawer here was the mistake: what this
% view is for is finding the cells two annotations disagree about, and the
% Multi-Grouping View already links the panels by brush, so brushing a
% suspicious blob in one panel lights up the same cells in every other. That
% is the question, answered directly. The camera is linked too, so a 3-D
% embedding stays comparable while it is rotated.
%
% What is given up is the per-panel legend: I_GSCATTER3 draws one SCATTER
% object coloured by group index, and a legend cannot be built from that. The
% toolbar's "Show group labels" button covers it by writing each type's name
% at its centroid, and hovering names the type under the cursor.

ntypes = cellfun(@(lbl) numel(unique(lbl)), labels);
titles = names(:) + " (" + pkg.i_plural(ntypes(:), 'type') + ")";
gui.i_multigroupview(sce, labels, titles, parentfig, ...
    'Cell Type Annotations');
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

indx = in_picktwo(names, parentfig, ['Select exactly two annotations. Rows ', ...
    'will be the first, columns the second.']);
if isempty(indx), return; end

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
%
% The adjusted Rand index is reported either way, and is the number to read
% when the vocabularies differ: it scores the two groupings on which cells
% they put together, never on what those groups are called.
ari = pkg.i_adjustedrandindex(a, b);
shared = intersect(la, lb);
if isempty(shared)
    msg = sprintf(['"%s" and "%s" share no label names, so a per-cell ', ...
        'agreement rate would read as 0%% however well the groupings line ', ...
        'up. Adjusted Rand index: %.3f (%s).'], ...
        names(indx(1)), names(indx(2)), ari, in_arireading(ari));
else
    msg = sprintf(['"%s" and "%s" give the same label to %.1f%% of cells ', ...
        '(%d of %d), over %d shared label name(s). Adjusted Rand index: ', ...
        '%.3f (%s).'], ...
        names(indx(1)), names(indx(2)), 100*mean(a == b), sum(a == b), ...
        numel(a), numel(shared), ari, in_arireading(ari));
end
gui.myHelpdlg(parentfig, msg);
end


function s = in_arireading(ari)
% A word for the number, because ARI has no intuitive scale: it is not a
% percentage, and the reference point that matters is 0 = chance, not 0 =
% no overlap. The bands are a reading aid, not a test.
if isnan(ari)
    s = 'not defined for fewer than two cells';
elseif ari >= 0.9
    s = '1 = same partition';
elseif ari >= 0.6
    s = 'largely the same partition, some types split or merged';
elseif ari >= 0.3
    s = 'partly overlapping partitions';
elseif ari > 0.05
    s = 'little more than chance agreement';
else
    s = '0 = chance agreement';
end
end
