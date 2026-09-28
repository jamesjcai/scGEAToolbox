function [p, W, zval] = e_ranksumrows(Z, grp, method)
%E_RANKSUMROWS Wilcoxon rank-sum test of every row, one group against the rest.
%   [P, W] = pkg.e_ranksumrows(Z, GRP) tests each row of Z (genes x cells)
%   for a difference between the cells of each group and all other cells.
%   GRP holds one positive integer label per cell (column of Z). P and W are
%   genes x numGroups: P(i,k) is the two-sided p-value for gene i, group k
%   against the rest, and W(i,k) is the rank sum of group k.
%
%   The numbers are RANKSUM's: P(i,k) equals ranksum(Z(i,GRP==k),
%   Z(i,GRP~=k)) and W(i,k) its STATS.RANKSUM, including the tie
%   correction, the continuity correction, the exact method RANKSUM picks
%   for small samples (fewest cells < 10 and all cells < 20), NaN for a row
%   with no variation, and NaNs dropped from each sample.
%
%   [P, W] = pkg.e_ranksumrows(Z, GRP, "approximate") uses the normal
%   approximation for every group, as ranksum(..., 'method', 'approximate').
%
%   [P, W, ZVAL] = ... also returns RANKSUM's STATS.ZVAL, the continuity-
%   corrected normal statistic, signed so that it is negative when group k
%   ranks low. RANKSUM reports no z for a group it tests exactly, and ZVAL
%   is NaN there too.
%
%   Why: calling RANKSUM once per gene re-sorts every cell of that gene,
%   although in a count matrix most of them are one tied block of zeros,
%   and a one-vs-rest marker search repeats that for every group. Here each
%   gene's nonzero values are ranked once and every group's rank sum is read
%   off those ranks, so the cost scales with the number of nonzeros, not with
%   genes x cells x groups, and a sparse Z is never made dense.
%
%   Example -- the two-sample test SC_DEG runs:
%       p = pkg.e_ranksumrows([X, Y], [ones(size(X,2),1); 2*ones(size(Y,2),1)]);
%       p = p(:, 1);
%
%   See also ranksum, tiedrank, sc_deg, pkg.e_findallmarkers.

arguments
    Z {mustBeNumeric}
    grp {mustBeVector, mustBePositive, mustBeInteger}
    method (1,1) string {mustBeMember(method, ["auto", "approximate"])} = "auto"
end

[G, N] = size(Z);
grp = double(grp(:));
if numel(grp) ~= N
    error("pkg:e_ranksumrows:SizeMismatch", ...
        "GRP has %d labels but Z has %d columns; give one label per cell.", ...
        numel(grp), N);
end
K = max(grp);
nk = accumarray(grp, 1, [K, 1]).';

p = nan(G, K);
W = nan(G, K);
zval = nan(G, K);

% Rows holding NaN go to RANKSUM itself: it drops NaNs per sample, which
% changes each row's sample sizes. Z(:) avoids building a dense logical of
% the whole matrix for a sparse Z.
nanrows = false(G, 1);
if any(isnan(nonzeros(Z)))
    [ri, ~] = find(isnan(Z));
    nanrows(ri) = true;
end

% Groups RANKSUM would test exactly (few cells) also go to RANKSUM.
small = min(nk, N - nk) < 10 & N < 20;
if method == "approximate"
    small(:) = false;
end

% Blocks of genes bound the memory of the nonzero lists (or, for a dense
% block, of the sorted copies).
Zt = Z.';                                   % cells x genes: column slices are fast
blockSize = max(1, floor(5e6/max(1, nnz(Z)/max(G, 1))));
for b0 = 1:blockSize:G
    rows = b0:min(G, b0 + blockSize - 1);
    [pb, Wb, zb] = i_block(Zt(:, rows), grp, nk, N, K);
    p(rows, :) = pb;
    W(rows, :) = Wb;
    zval(rows, :) = zb;
end

% The fallbacks, row by row and group by group.
exactGroups = find(small);
for i = 1:G
    if nanrows(i)
        ks = 1:K;
    elseif ~isempty(exactGroups)
        ks = exactGroups;
    else
        continue;
    end
    z = full(Z(i, :));
    for k = ks
        if nk(k) == 0 || nk(k) == N, continue; end
        if method == "approximate"
            [p(i, k), ~, st] = ranksum(z(grp == k), z(grp ~= k), 'method', 'approximate');
        else
            [p(i, k), ~, st] = ranksum(z(grp == k), z(grp ~= k));
        end
        W(i, k) = st.ranksum;
        if isfield(st, 'zval') && ~isempty(st.zval)
            zval(i, k) = st.zval;
        else
            zval(i, k) = NaN;       % the exact test reports no z
        end
    end
end
end


function [p, W, zval] = i_block(Zt, grp, nk, N, K)
% Normal-approximation p-values and rank sums for one block of genes.
% ZT is cells x genes.
nG = size(Zt, 2);
if nnz(Zt) > 0.25*numel(Zt)
    % Mostly nonzero (Pearson residuals, z-scores): sorting each column
    % beats listing and sorting (gene, value) pairs.
    [W, tieadj] = i_denseranks(full(double(Zt)), grp, N, K);
else
    [W, tieadj] = i_sparseranks(Zt, grp, nk, N, K);
end

% RANKSUM's normal approximation, from the smaller sample's rank sum.
nx = nk;                                     % the group
ny = N - nk;                                 % the rest
same = nx <= ny;
ns = min(nx, ny);
w = W;
w(:, ~same) = N*(N + 1)/2 - W(:, ~same);     % rank sum of the rest
wmean = ns*(N + 1)/2;
tiescor = 2*tieadj/(N*(N - 1));
wvar = nx.*ny.*((N + 1) - tiescor)/12;
wc = w - wmean;
z = (wc - 0.5*sign(wc))./sqrt(wvar);
p = 2*normcdf(-abs(z));
p(:, nk == 0 | nk == N) = NaN;
% RANKSUM's STATS.ZVAL is z for the first sample: flipped when the group was
% the larger sample and z was computed from the rest's rank sum.
zval = z;
zval(:, ~same) = -z(:, ~same);
zval(:, nk == 0 | nk == N) = NaN;
end


function [W, tieadj] = i_denseranks(D, grp, N, K)
% Rank sums per group and TIEDRANK's tie adjustment for every column of D
% (cells x genes), ties averaged.
nG = size(D, 2);
[vs, o] = sort(D, 1);
pos = repmat((1:N).', 1, nG);
newtie = [true(1, nG); diff(vs, 1, 1) ~= 0];
lasttie = [diff(vs, 1, 1) ~= 0; true(1, nG)];
first = pos;
first(~newtie) = 0;
first = cummax(first, 1);
last = pos;
last(~lasttie) = Inf;
last = flipud(cummin(flipud(last), 1));
R = zeros(N, nG);
R(o + (0:nG-1)*N) = (first + last)/2;
% Ranks are multiples of 1/2, so these sums are exact in any order.
W = full(R.' * sparse(1:N, grp, 1, N, K));
tsize = last - first + 1;
tieadj = (sum((tsize.^3 - tsize).*newtie, 1)/2).';
end


function [W, tieadj] = i_sparseranks(Zt, grp, nk, N, K)
% Rank sums per group and the tie adjustment from the stored entries only.
nG = size(Zt, 2);
[j, i, v] = find(Zt);                        % j cell, i gene within block
v = full(double(v));

% Sort each gene's nonzero values. FIND returns them by gene (column of
% ZT), so a stable sort on value within gene is a sort on [i, v].
[~, o] = sortrows([i, v]);
i = i(o);
j = j(o);
v = v(o);

nnzg = accumarray(i, 1, [nG, 1]);
nneg = accumarray(i, v < 0, [nG, 1]);
n0 = N - nnzg;                               % the tied block of zeros

% Position of each nonzero among its gene's nonzeros, averaged over ties.
first = cumsum([0; nnzg(1:end-1)]);
pos = (1:numel(i)).' - first(i);
newtie = [true; diff(i) ~= 0 | diff(v) ~= 0];
tid = cumsum(newtie);
tsize = accumarray(tid, 1);
tfirst = pos(newtie);
avgpos = tfirst(tid) + (tsize(tid) - 1)/2;

% Global rank among all N cells: negatives rank first, then the zeros,
% then the positives.
r = avgpos + (v > 0).*n0(i);
zeroRank = nneg + (n0 + 1)/2;

% Rank sum of every group: its nonzeros plus its zeros at the zero rank.
Wnz = accumarray([i, grp(j)], r, [nG, K]);
cnt = accumarray([i, grp(j)], 1, [nG, K]);
W = Wnz + (nk - cnt).*zeroRank;

% Tie adjustment, as TIEDRANK returns it: sum(t^3 - t)/2 over tie groups,
% the zeros being one group.
tieg = i(newtie);
tieadj = accumarray(tieg, (tsize.^3 - tsize)/2, [nG, 1]) + (n0.^3 - n0)/2;
end
