function hFig = i_alluvialview(a, b, nameA, nameB, parentfig, figname)
%I_ALLUVIALVIEW Sankey/alluvial diagram between two groupings of the same cells.
%
%   hFig = gui.i_alluvialview(a, b, nameA, nameB)
%   hFig = gui.i_alluvialview(a, b, nameA, nameB, parentfig, figname)
%
%   A and B are per-cell label vectors of equal length - two cell type
%   annotations, a clustering and an annotation, before and after a merge.
%   Each is drawn as a column of blocks, one per type, sized by cell count,
%   and a ribbon carries each group of cells from its label in A to its label
%   in B.
%
%   This answers the question the cross-tabulation answers, in the form the
%   question is usually asked in: which type got split, which two got merged,
%   and where did a type that vanished go. A table of the same numbers makes
%   the reader do that reconstruction in their head.
%
%   Ribbons are coloured by their source, so a split shows as one colour
%   fanning out and a merge as several colours converging. Right-hand blocks
%   are deliberately neutral: colouring them from the same palette would put
%   the same colour on the left and right of a diagram whose whole point is
%   that the two vocabularies need not line up.
%
%   Click a ribbon or a block for its counts.
%
%   Returns the figure handle, or [] if the two vectors cannot be compared.
%
%   See also PKG.I_ALLUVIALLAYOUT, PKG.I_ADJUSTEDRANDINDEX,
%   GUI.CALLBACK_COMPARECELLTYPEANNOTATIONS.

if nargin < 6, figname = 'Annotation Flow'; end
if nargin < 5, parentfig = []; end
if nargin < 4 || isempty(nameB), nameB = "B"; end
if nargin < 3 || isempty(nameA), nameA = "A"; end

hFig = [];
a = string(a(:));
b = string(b(:));
if numel(a) ~= numel(b)
    error('gui:i_alluvialview:sizeMismatch', ...
        ['The two label vectors have %d and %d elements. Both must label ', ...
        'the same cells.'], numel(a), numel(b));
end
if isempty(a), return; end

[ga, la] = findgroups(a);
[gb, lb] = findgroups(b);
M = accumarray([ga, gb], 1, [numel(la), numel(lb)]);
layout = pkg.i_alluviallayout(M);

nA = numel(la);
nB = numel(lb);
cmap = pkg.i_mycolorlines(nA);

hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;
hFig.Name = char(figname);
hFig.NumberTitle = 'off';

% GUI.MYFIGURE already made an axes, and its toolbar's save and camera
% buttons were handed that one. Drawing into a second axes left the first
% underneath with its default 0-1 rulers showing through.
ax = hx.AxHandle;
cla(ax, 'reset');
hold(ax, 'on');

xLeft = 0;
xRight = 1;
barWidth = 0.055;

% Ribbons first, so the blocks sit on top of where they meet.
rb = layout.ribbon;
for k = 1:numel(rb.count)
    [xs, yTop, yBot] = in_ribbon(xLeft + barWidth, xRight - barWidth, ...
        rb.lTop(k), rb.lBot(k), rb.rTop(k), rb.rBot(k));
    p = patch(ax, [xs, fliplr(xs)], [yTop, fliplr(yBot)], ...
        cmap(rb.i(k), :), 'EdgeColor', 'none', 'FaceAlpha', 0.45);
    p.UserData = sprintf('%s -> %s: %s (%.1f%% of %s)', ...
        la(rb.i(k)), lb(rb.j(k)), pkg.i_plural(rb.count(k), 'cell'), ...
        100*rb.count(k)/sum(M(rb.i(k), :)), la(rb.i(k)));
end

% Blocks and their labels.
for i = 1:nA
    in_block(ax, xLeft, barWidth, layout.leftTop(i), layout.leftBot(i), ...
        cmap(i, :), sprintf('%s: %s', la(i), ...
        pkg.i_plural(sum(M(i, :)), 'cell')));
end
for j = 1:nB
    in_block(ax, xRight - barWidth, barWidth, layout.rightTop(j), ...
        layout.rightBot(j), [0.6 0.6 0.6], ...
        sprintf('%s: %s', lb(j), pkg.i_plural(sum(M(:, j)), 'cell')));
end

nshown = in_labels(ax, la, layout.leftTop, layout.leftBot, xLeft - 0.015, ...
    'right', layout.gap);
nshown = nshown + in_labels(ax, lb, layout.rightTop, layout.rightBot, ...
    xRight + 0.015, 'left', layout.gap);

% Column headings, at a y clear of the tallest block.
text(ax, xLeft + barWidth/2, 1.045, nameA, 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'Interpreter', 'none');
text(ax, xRight - barWidth/2, 1.045, nameB, 'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', 'Interpreter', 'none');

ari = pkg.i_adjustedrandindex(a, b);
% Two independent groupings land just either side of zero, and "-0.000"
% reads as a real negative rather than as chance.
if abs(ari) < 0.0005, ari = 0; end
ttl = sprintf('%d cells, %d -> %d types, ARI %.3f', numel(a), nA, nB, ari);
if nshown < nA + nB
    % Say it rather than leave the reader counting blocks against labels.
    ttl = sprintf('%s   (%d thin block(s) unlabelled - click for counts)', ...
        ttl, nA + nB - nshown);
end
title(ax, ttl, 'Interpreter', 'none');

axis(ax, 'off');
% Room either side for the labels, which are drawn outside the columns.
xlim(ax, [-0.32, 1.32]);
ylim(ax, [-0.03, 1.09]);

dt = datacursormode(hFig);
dt.UpdateFcn = @in_datatip;
hx.show(parentfig);
end


function [xs, yTop, yBot] = in_ribbon(x1, x2, lTop, lBot, rTop, rBot)
% A cosine ease from one column to the other: flat where it meets each block,
% steepest in the middle. The two edges share it, so the ribbon keeps its
% thickness measured vertically and a thin flow stays visible along its whole
% length.
npt = 40;
t = linspace(0, 1, npt);
ease = 0.5*(1 - cos(pi*t));
xs = x1 + (x2 - x1)*t;
yTop = lTop + (rTop - lTop)*ease;
yBot = lBot + (rBot - lBot)*ease;
end


function in_block(ax, x, w, top, bot, color, info)
p = patch(ax, [x, x+w, x+w, x], [top, top, bot, bot], color, ...
    'EdgeColor', 'none');
p.UserData = info;
end


function nshown = in_labels(ax, names, top, bot, x, align, gap)
% A label per block, skipped where the block is thinner than the gap between
% blocks - below that the text overlaps its neighbours and the column becomes
% unreadable. The title reports how many were skipped.
nshown = 0;
for k = 1:numel(names)
    if top(k) - bot(k) < gap, continue; end
    text(ax, x, (top(k) + bot(k))/2, names(k), ...
        'HorizontalAlignment', align, 'VerticalAlignment', 'middle', ...
        'FontSize', 8, 'Interpreter', 'none');
    nshown = nshown + 1;
end
end


function txt = in_datatip(target, ~)
% Every patch carries its own description, so one function serves ribbons and
% blocks alike.
txt = '';
if isprop(target, 'UserData') && ~isempty(target.UserData)
    txt = gui.i_escapeunderscore(string(target.UserData));
end
end
