function i_qcviolin(X, genelist, parentfig, groups, groupname, grouporder)
%I_QCVIOLIN Violin plots of the three cell QC metrics.
%
%   gui.i_qcviolin(X, genelist, parentfig) plots genes detected, library
%   size and mitochondrial percentage over all cells, side by side.
%
%   gui.i_qcviolin(X, genelist, parentfig, groups, groupname) plots one
%   violin per group instead, the three metrics stacked so the groups have
%   the width. GROUPS has one entry per cell; GROUPNAME labels the x-axis.
%
%   gui.i_qcviolin(..., grouporder) draws the groups in the order of
%   GROUPORDER, the level names; without it they are sorted.

if nargin < 3, parentfig = []; end
if nargin < 4, groups = []; end
if nargin < 5, groupname = ''; end
if nargin < 6, grouporder = []; end

i = startsWith(genelist, 'mt-', 'IgnoreCase', true);
nftr = full(sum(X > 0, 1));
lbsz = full(sum(X, 1));
lbsz_mt = full(sum(X(i, :), 1));
cj = 100 * (lbsz_mt ./ lbsz);

metrics = {nftr, lbsz, cj};
titles = {sprintf('nFeature\\_RNA\n(# of genes)'), ...
    sprintf('nCount\\_RNA\n(# of reads)'), ...
    sprintf('percent.mt\n(mitochondrial content)')};

hx = gui.myFigure(parentfig, true);

if isempty(groups)
    for k = 1:3
        subplot(1, 3, k)
        gui.i_violinplot_base(metrics{k}, [], 'showdata', false);
        title(titles{k});
        box on;
    end
else
    % Three stacked rows and long rotated group names do not fit the
    % default figure: the axes were squeezed to a line. Make it taller,
    % pack the tiles, and label the groups under the bottom row only.
    groupedFigureHeight = 760;
    widthPerGroup = 55;
    widthMargin = 200;
    fx = hx.FigHandle;
    numGroups = numel(unique(string(groups)));
    fx.Position(3) = max(fx.Position(3), ...
        min(widthPerGroup*numGroups + widthMargin, 1400));
    fx.Position(4) = groupedFigureHeight;
    t = tiledlayout(fx, 3, 1, TileSpacing="compact", Padding="compact");
    % Underscores would be read as TeX subscripts in the tick labels.
    cats = strrep(string(groups(:)), '_', ' ');
    orderargs = {};
    if ~isempty(grouporder)
        orderargs = {'GroupOrder', cellstr(strrep(string(grouporder(:)), '_', ' '))};
    end
    for k = 1:3
        ax = nexttile(t);
        gui.i_violinplot_base(metrics{k}(:), cats, 'showdata', false, orderargs{:});
        ylabel(ax, titles{k});
        box(ax, "on");
        if k < 3
            ax.XTickLabel = [];
        else
            xtickangle(ax, -45);
        end
    end
    if strlength(string(groupname)) > 0
        xlabel(t, strrep(string(groupname), '_', ' '));
    end
end

hx.show(parentfig);

end
