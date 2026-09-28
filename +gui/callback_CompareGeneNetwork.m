function callback_CompareGeneNetwork(src, ~)
% CALLBACK_COMPAREGENENETWORK  Build a GRN for each of two cell groups and compare.
%
% One flow for Network > Build && Compare Two GRNs...:
%
%   1. Cells   all cells, or the cells in chosen groups of one grouping
%   2. Groups  a grouping, then the group(s) on each side of the comparison;
%              skipped when the grouping has exactly two groups
%   3. Genes   pick from a list, paste names, or all genes
%   4. Method  selected genes: any method in net.grnmethods, then the
%              transform, and the two networks are drawn side by side with
%              sc_grnview2. All genes: gui.callback_CompareGRNAllGenes runs
%              the scTenifoldNet comparison.
%
% See also gui.callback_CompareGRNAllGenes, gui.callback_BuildGeneNetwork,
%   net.grnmethods, sc_grnview2.

[FigureHandle, sce] = gui.gui_getfigsce(src);

% --- Cells
cellIdx = gui.i_pickgrncells(sce, FigureHandle, sprintf(['Two networks ', ...
    'are best compared within one cell population, e.g. one cell type ', ...
    'in two conditions. Compare within all %d cells, or only within ', ...
    'cells in groups you select?'], sce.NumCells));
if isempty(cellIdx), return; end

% --- Groups
[i1, i2, name1, name2] = i_picktwogroups(sce, cellIdx, FigureHandle);
if isempty(i1), return; end

% --- Genes
[glist, useAll] = gui.i_pickgrngenes(sce.g, FigureHandle);
if useAll
    gui.callback_CompareGRNAllGenes(src, i1, i2, name1, name2);
    return;
end
if isempty(glist), return; end

[y, i] = ismember(upper(glist), upper(sce.g));
if ~all(y), error('Selected gene(s) not in the gene list of data.'); end
fprintf("%s\n", glist)

% The method list, its preferred transform and its pre-run warnings all
% live in net.grnmethods. The transform used to be preselected here as
% item 5, which meant Pearson residuals until "(d): Shifted CLR" was
% inserted into gui.i_transformx's list and made it DESeq.
method = gui.i_selectgrnmethod(FigureHandle, sprintf(['How should the ', ...
    'links between the %d genes be scored? The same method builds both ', ...
    'networks.'], numel(glist)));
if isempty(method), return; end

% Transform only the cells being compared, together, so both networks
% see the same scale.
inEither = i1 | i2;
[Xt] = gui.i_transformx(sce.X(:, inEither), true, method.Transform, FigureHandle);
if isempty(Xt), return; end
if ~gui.i_confirmgrnrun(method, Xt(i, :), 2, FigureHandle), return; end

x1 = Xt(i, i1(inEither));
x2 = Xt(i, i2(inEither));

fw = gui.myWaitbar(FigureHandle);
try
    A1 = sc_grn(x1, method.Key);
    A2 = sc_grn(x2, method.Key);
catch ME
    gui.myWaitbar(FigureHandle, fw, true);
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
    return;
end
gui.myWaitbar(FigureHandle, fw);

drawnow;
stitle = sprintf('%s vs. %s', name1, name2);
try
    sc_grnview2(A1, A2, glist, stitle, FigureHandle);
catch ME
    gui.myErrordlg(FigureHandle, ME.message, ME.identifier);
end
end


function [i1, i2, name1, name2] = i_picktwogroups(sce, cellIdx, parentfig)
% Logical indices (one per cell of SCE) of the two groups to compare,
% both inside CELLIDX, and their display names; all empty if cancelled.
i1 = [];
i2 = [];
name1 = "";
name2 = "";

[names, values] = gui.i_cellgroupings(sce, cellIdx);
if isempty(names)
    gui.myWarndlg(parentfig, ['Every grouping puts these cells in a ', ...
        'single group, so there are no two groups to compare. Choose ', ...
        'more cells, or add a grouping such as a batch or condition.']);
    return;
end
if isscalar(names)
    k = 1;
else
    [k, ok] = gui.myListdlg(parentfig, names, 'Compare Networks By', 1, ...
        false, true, [300, 300], ['Select the grouping that defines the ', ...
        'two groups to compare.']);
    if ~ok || isempty(k), return; end
end

v = values{k};
v(~cellIdx) = missing;
[groups, labels] = gui.i_grouplabels(v(cellIdx));

if numel(groups) == 2
    % Nothing to choose; keep the groups in natural order.
    [~, order] = natsort(groups);
    sel1 = order(1);
    sel2 = order(2);
else
    [sel1, ok] = gui.myListdlg(parentfig, labels, 'Group 1', [], true, ...
        true, [340, 400], sprintf(['Select the %s group(s) that form ', ...
        'the first side of the comparison.'], names(k)));
    if ~ok || isempty(sel1), return; end
    rest = setdiff(1:numel(groups), sel1);
    if isempty(rest)
        gui.myWarndlg(parentfig, ['Every group is on the first side, so ', ...
            'none is left to compare with. Leave at least one out.']);
        return;
    end
    [pick, ok] = gui.myListdlg(parentfig, labels(rest), 'Group 2', [], ...
        true, true, [340, 400], sprintf(['Select the %s group(s) to ', ...
        'compare with %s.'], names(k), strjoin(groups(sel1), " + ")));
    if ~ok || isempty(pick), return; end
    sel2 = rest(pick);
end

i1 = ismember(v, groups(sel1));
i2 = ismember(v, groups(sel2));
name1 = strjoin(groups(sel1), " + ");
name2 = strjoin(groups(sel2), " + ");
end
