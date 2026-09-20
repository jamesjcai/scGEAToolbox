function layout = i_alluviallayout(M, gapEach)
%I_ALLUVIALLAYOUT Block and ribbon geometry for a two-column alluvial diagram.
%
%   layout = pkg.i_alluviallayout(M)
%   layout = pkg.i_alluviallayout(M, gapEach)
%
%   M is a contingency matrix: M(i,j) items are in left group I and right
%   group J. GAPEACH is the vertical gap between two neighbouring blocks, as a
%   fraction of the axis height (default 0.01, shrunk automatically when there
%   are too many groups for the gaps to fit).
%
%   Returns a struct whose fields are all indexed by the ORIGINAL group
%   numbers, i.e. M's own row and column indices, not the display order:
%
%     leftTop, leftBot     nA x 1, y of each left block's top and bottom
%     rightTop, rightBot   nB x 1, likewise on the right
%     ordLeft, ordRight    display order, top to bottom
%     ribbon               struct of arrays, one entry per nonzero M(i,j):
%                            .i .j .count and the four y values .lTop .lBot
%                            .rTop .rBot where the ribbon meets each column
%
%   y runs 0 at the bottom to 1 at the top. Both columns share one scale, so a
%   ribbon is the same thickness at both ends and thickness is comparable
%   across the diagram; the shorter column is centred rather than stretched.
%
%   Two ordering decisions do the work that makes such a diagram readable:
%
%   1. Left blocks are ordered by size, largest at the top. Right blocks are
%      then ordered by the barycentre of their incoming flow - the average
%      height of where their cells come from - which is the standard one-pass
%      crossing reduction, and puts a right group directly across from
%      whatever fed it.
%   2. Within a block, the ribbons are stacked in the order of the blocks they
%      connect to. Without this, ribbons cross each other inside a block even
%      when the blocks themselves are well ordered.
%
%   Geometry only: it draws nothing and needs no figure, which is what makes
%   it testable. See GUI.I_ALLUVIALVIEW for the drawing.
%
%   See also GUI.I_ALLUVIALVIEW, GUI.CALLBACK_COMPARECELLTYPEANNOTATIONS.

if nargin < 2 || isempty(gapEach), gapEach = 0.01; end

M = full(double(M));
[nA, nB] = size(M);
total = sum(M(:));
if total <= 0
    error('pkg:i_alluviallayout:emptyMatrix', ...
        'M holds no items. Pass a contingency matrix with a positive total.');
end

% One gap size for both columns, shrunk if the column with more blocks would
% otherwise spend more than a quarter of the height on gaps.
maxGaps = max(nA, nB) - 1;
if maxGaps > 0
    gapEach = min(gapEach, 0.25/maxGaps);
else
    gapEach = 0;
end

rowTot = sum(M, 2);
colTot = sum(M, 1)';

% The taller column fills the axis; the shared scale follows from it.
scale = (1 - gapEach*maxGaps)/total;

[~, ordLeft] = sort(rowTot, 'descend');
[leftTop, leftBot] = in_stack(rowTot, ordLeft, scale, gapEach);

% Barycentre of each right block's sources, in the y the left column just
% got. A right block with no incoming flow cannot have one, and is sent to
% the bottom rather than left as NaN, where SORT would scatter it.
leftMid = (leftTop + leftBot)/2;
bary = -inf(nB, 1);
for j = 1:nB
    w = M(:, j);
    if sum(w) > 0
        bary(j) = sum(w.*leftMid)/sum(w);
    end
end
[~, ordRight] = sort(bary, 'descend');
[rightTop, rightBot] = in_stack(colTot, ordRight, scale, gapEach);

% Ribbons. Each block is filled from its top downwards, taking its
% connections in the display order of the blocks at the other end, so no two
% ribbons cross within a block.
rightMid = (rightTop + rightBot)/2;
rankLeft = zeros(nA, 1); rankLeft(ordLeft) = 1:nA;
rankRight = zeros(nB, 1); rankRight(ordRight) = 1:nB;

% Explicit column shapes throughout. FIND on a matrix with a single row or
% column returns rows, and indexing a vector takes the orientation of the
% vector being indexed, so a one-type annotation produced a row on one side
% and a column on the other and the two would not concatenate.
[ii, jj] = find(M > 0);
ii = ii(:);
jj = jj(:);
counts = reshape(M(sub2ind([nA, nB], ii, jj)), [], 1);
nrib = numel(ii);
lTop = zeros(nrib, 1); lBot = zeros(nrib, 1);
rTop = zeros(nrib, 1); rBot = zeros(nrib, 1);

rankOfLeft = reshape(rankLeft(ii), [], 1);
rankOfRight = reshape(rankRight(jj), [], 1);

offset = leftTop;
[~, bySource] = sortrows([rankOfLeft, rankOfRight]);
for k = reshape(bySource, 1, [])
    lTop(k) = offset(ii(k));
    lBot(k) = lTop(k) - counts(k)*scale;
    offset(ii(k)) = lBot(k);
end

offset = rightTop;
[~, byTarget] = sortrows([rankOfRight, rankOfLeft]);
for k = reshape(byTarget, 1, [])
    rTop(k) = offset(jj(k));
    rBot(k) = rTop(k) - counts(k)*scale;
    offset(jj(k)) = rBot(k);
end

layout = struct( ...
    'leftTop', leftTop, 'leftBot', leftBot, ...
    'rightTop', rightTop, 'rightBot', rightBot, ...
    'leftMid', leftMid, 'rightMid', rightMid, ...
    'ordLeft', ordLeft, 'ordRight', ordRight, ...
    'scale', scale, 'gap', gapEach, 'total', total, ...
    'ribbon', struct('i', ii, 'j', jj, 'count', counts, ...
                     'lTop', lTop, 'lBot', lBot, ...
                     'rTop', rTop, 'rBot', rBot));
end


function [top, bot] = in_stack(tot, ord, scale, gapEach)
% Stack one column's blocks top to bottom in display order ORD, then centre
% the result vertically. Returned in original index order.
n = numel(tot);
top = zeros(n, 1);
bot = zeros(n, 1);
y = 1;
for k = 1:n
    g = ord(k);
    top(g) = y;
    bot(g) = y - tot(g)*scale;
    y = bot(g) - gapEach;
end

used = 1 - (y + gapEach);
shift = (1 - used)/2;
top = top - shift;
bot = bot - shift;
end
