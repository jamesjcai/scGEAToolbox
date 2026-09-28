function callback_CompareCellTypeAnnotations(src, ~)
%CALLBACK_COMPARECELLTYPEANNOTATIONS Compare cell type annotations in a Sankey diagram.
%
%   Every annotation channel stashes the labels it replaces as a numbered
%   'old_cell_type_N' cell attribute (PKG.I_STASHCELLTYPEHISTORY), so a dataset
%   annotated more than once carries a history. Keeping it is only half the
%   job: this is how it gets looked at.
%
%   The Sankey diagram between two annotations is what opens: which types
%   got split or merged is the question this is asked most often, and the
%   flow answers it at a glance. The other views are buttons on its toolbar,
%   because the useful question differs:
%     - side-by-side embedding plots, linked by brushing, answer "where do
%       these disagree", across every annotation, not just the two shown;
%     - the per-cell table answers "what did each method call THIS cell";
%     - the cross-tabulation gives the exact counts behind the ribbons, for
%       the same two annotations the diagram shows.
%
%   See also PKG.I_CELLTYPEHISTORY, PKG.I_STASHCELLTYPEHISTORY,
%   GUI.I_UPDATEANNOTATEMENU, GUI.I_ALLUVIALVIEW.

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

% With exactly two there is only one pair, so asking would be a formality.
if numel(names) == 2
    pick = [1, 2];
else
    pick = in_picktwo(names, parentfig, ['Select exactly two annotations. ', ...
        'The flow runs from the first to the second.']);
    if isempty(pick), return; end
end

% Views opened from the toolbar are parented to the Sankey window, so they
% open beside the diagram they were asked from, not over the main window.
% The panel view is drawn into the Sankey window itself, with a Back button.
% The tables are uifigures and cannot share it, so they open in its place on
% screen instead, with a Back button to it.
buttons = struct( ...
    'Icon', {'brush.gif', ...
    'data_table_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', ...
    'view-grid.jpg'}, ...
    'Tooltip', {'Side-by-side plots of all annotations (linked by brushing)', ...
    'Table with one row per cell, all annotations', ...
    'Cross-tabulate these two annotations (exact counts)'}, ...
    'Callback', { ...
    @(src, ~) in_plotpanels(sce, names, labels, ancestor(src, 'figure')), ...
    @(src, ~) in_showtable(sce, names, labels, ancestor(src, 'figure')), ...
    @(src, ~) in_crosstabshown(ancestor(src, 'figure'))});

gui.i_alluvialview(labels{pick(1)}, labels{pick(2)}, names(pick(1)), ...
    names(pick(2)), parentfig, 'Cell Type Annotation Flow', buttons);
end


function pick = in_picktwo(names, parentfig, prompt)
% Which pair the Sankey diagram shows; the cross-tabulation button reuses it.
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
gui.myFigure.drawInto(parentfig, @() gui.i_multigroupview(sce, labels, ...
    titles, parentfig, 'Cell Type Annotations'));
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

gui.i_openinstead(parentfig, @() gui.TableViewerApp(t, parentfig, 'CellTypeAnnotations'));
end


function in_crosstabshown(hFig)
% Cross-tabulate the pair in the direction the Sankey diagram currently
% draws it - its swap button may have reversed it since it opened - so the
% rows are always the left column and the columns the right one.
P = getappdata(hFig, 'AlluvialPair');
in_crosstab([P.NameA, P.NameB], {P.A, P.B}, hFig);
end


function in_crosstab(names, labels, parentfig)
% Cross-tabulate two annotations: rows the first, columns the second, counts
% inside. NAMES and LABELS hold just the pair the Sankey diagram shows.

indx = [1, 2];
a = labels{1};
b = labels{2};
[ga, la] = findgroups(a);
[gb, lb] = findgroups(b);
M = accumarray([ga(:), gb(:)], 1, [numel(la), numel(lb)]);

colnames = matlab.lang.makeUniqueStrings(matlab.lang.makeValidName(lb));
t = array2table(M, 'VariableNames', colnames);
t = addvars(t, la(:), 'Before', 1, 'NewVariableNames', {'RowLabel'});
gui.i_openinstead(parentfig, @() gui.TableViewerApp(t, parentfig, 'CellTypeCrosstab'));

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
