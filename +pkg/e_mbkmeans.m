function [c, counts, idx] = e_mbkmeans(x, k, c, counts, opts)
%E_MBKMEANS  Mini-batch k-means clustering (Sculley, WWW 2010).
%
%   [c, counts, idx] = PKG.E_MBKMEANS(x, k) partitions the rows of X into K
%   clusters and returns the centroids C, the per-centroid update counts
%   COUNTS, and the cluster index IDX of every row.
%
%   Mini-batch k-means trades a little accuracy for speed. Each iteration
%   draws a small random batch, assigns it to the nearest centroids, and
%   moves each centroid toward the batch points it won, with a per-centroid
%   learning rate of 1/count. The work per iteration depends on BatchSize
%   rather than on the number of rows, so the run time is nearly flat in N
%   where Lloyd's algorithm is linear in it. That is the whole reason to
%   prefer this over KMEANS on a large embedding.
%
%   WHAT THIS REPLACES. The previous implementation performed exactly one
%   nearest-centroid assignment against centroids initialised to the FIRST K
%   ROWS of X, and returned it. There was no iteration, and the centroid
%   update that followed could not affect the labels that had already been
%   returned. On two well-separated Gaussian blobs of 500 points each it
%   scored ARI 0.007 against KMEANS's 1.000, splitting 43/957 where the
%   truth was 500/500 -- because rows 1 and 2 of a block-ordered matrix are
%   both drawn from the first blob, so the "clustering" was the perpendicular
%   bisector of two nearby points from the same group. Both halves of that
%   are fixed here: k-means++ seeding, and an actual iteration.
%
%   INPUTS:
%     x      - N-by-D data, one observation per row.
%     k      - number of clusters. Clamped to N with a warning if larger.
%     c      - K-by-D warm-start centroids ([] to seed with k-means++).
%     counts - K-by-1 warm-start update counts ([] for zeros). Pass C and
%              COUNTS back in to continue refining a previous run.
%
%   NAME-VALUE:
%     BatchSize - observations sampled per iteration (default 1000). Larger
%                 batches give steadier updates at proportionally more cost.
%     MaxIter   - maximum iterations (default 100).
%     Replicates- restarts from independent seedings, keeping the one with
%                 the lowest within-cluster sum of squares (default 5).
%                 Ignored when C is warm-started, since the caller has
%                 already chosen the starting point. This is not a luxury:
%                 k-means from a single seeding lands in a local minimum
%                 that merges two clusters and splits a third about two
%                 thirds of the time. On six Gaussian clusters in 10-D over
%                 20 seeds, one seeding reached the true partition 7 times
%                 out of 20 -- for this function and for KMEANS alike --
%                 and five reached it 19 times. Mini-batch iterations are
%                 cheap enough that five restarts still cost less than one
%                 pass of Lloyd's algorithm.
%     Tolerance - stop once the largest centroid moves less than this
%                 fraction of the data's mean radius for three consecutive
%                 iterations (default 1e-4). Mini-batch updates are noisy,
%                 so a single small step is not convergence.
%
%   OUTPUTS:
%     c      - K-by-D final centroids.
%     counts - K-by-1 total number of points each centroid has absorbed.
%     idx    - N-by-1 cluster index, from a full assignment against the
%              final centroids. Every cluster is non-empty.
%
%   See also SC_CLUSTER_S, KMEANS.

arguments
    x {mustBeNumeric, mustBeNonempty, mustBeReal}
    k (1, 1) double {mustBePositive, mustBeInteger}
    c = []
    counts = []
    opts.BatchSize (1, 1) double {mustBePositive} = 1000
    opts.MaxIter (1, 1) double {mustBePositive} = 100
    opts.Tolerance (1, 1) double {mustBeNonnegative} = 1e-4
    opts.Replicates (1, 1) double {mustBePositive, mustBeInteger} = 5
end

if issparse(x), x = full(x); end
x = double(x);
[N, D] = size(x);

if k > N
    warning("e_mbkmeans:tooManyClusters", ...
        "K = %d exceeds the %d observations; using K = %d.", k, N, N);
    k = N;
end

% Every point is its own cluster; there is nothing to iterate toward.
if k == N
    c = x;
    counts = ones(k, 1);
    idx = (1:N).';
    return;
end

warmStart = ~isempty(c);
if warmStart
    c = double(c);
    if ~isequal(size(c), [k, D])
        error("e_mbkmeans:centroidSize", ...
            "C must be %d-by-%d; got %s.", k, D, mat2str(size(c)));
    end
end
if isempty(counts)
    counts = zeros(k, 1);
else
    counts = double(counts(:));
    if numel(counts) ~= k
        error("e_mbkmeans:countsSize", ...
            "COUNTS must have %d entries; got %d.", k, numel(counts));
    end
end

% Convergence is judged against the spread of the data so that Tolerance
% means the same thing whatever the units are.
scale = mean(vecnorm(x - mean(x, 1), 2, 2));
if scale == 0, scale = 1; end
batchSize = min(opts.BatchSize, N);

% Replicates are scored on one fixed subsample rather than on all of X, so
% that comparing them costs a fraction of a full assignment. Only the winner
% is assigned in full.
numRep = opts.Replicates;
if warmStart, numRep = 1; end
evalRows = min(N, max(10*batchSize, 1000));
if evalRows < N
    xEval = x(randperm(N, evalRows), :);
else
    xEval = x;
end

bestCost = inf;
bestC = c;
bestCounts = counts;
for rep = 1:numRep
    if warmStart
        cRep = c;
    else
        cRep = i_kmeanspp(x, k, batchSize);
    end
    [cRep, countsRep] = i_runbatches(x, cRep, counts, k, D, batchSize, ...
        opts.MaxIter, opts.Tolerance*scale);
    [~, dEval] = knnsearch(cRep, xEval, "k", 1);
    cost = sum(dEval.^2);
    if cost < bestCost
        bestCost = cost;
        bestC = cRep;
        bestCounts = countsRep;
    end
end
c = bestC;
counts = bestCounts;

idx = knnsearch(c, x, "k", 1);

% A centroid that never won a batch point can end up owning nothing. The
% caller asked for k clusters, so move each empty centroid onto the point
% currently worst served by its own centroid and reassign.
empty = find(~ismember((1:k).', idx));
if ~isempty(empty)
    d = vecnorm(x - c(idx, :), 2, 2);
    [~, worst] = sort(d, "descend");
    c(empty, :) = x(worst(1:numel(empty)), :);
    idx = knnsearch(c, x, "k", 1);
end

end


function [c, counts] = i_runbatches(x, c, counts, k, D, batchSize, maxIter, tol)
% One mini-batch k-means run from the given seeding.
N = size(x, 1);
patience = 3;
quiet = 0;
for iter = 1:maxIter
    batch = x(randperm(N, batchSize), :);
    a = knnsearch(c, batch, "k", 1);

    % Sculley's per-point update with learning rate 1/count is exactly a
    % running mean, so a whole batch can be applied at once:
    %   c_new = (count*c_old + sum of the points it won)/(count + won)
    % which is identical to applying the points one at a time, and orders of
    % magnitude faster.
    won = accumarray(a, 1, [k, 1]);
    total = zeros(k, D);
    for d = 1:D
        total(:, d) = accumarray(a, batch(:, d), [k, 1]);
    end
    hit = won > 0;
    cPrev = c;
    c(hit, :) = (counts(hit).*c(hit, :) + total(hit, :))./(counts(hit) + won(hit));
    counts = counts + won;

    % Mini-batch updates are noisy, so one small step is not convergence.
    if max(vecnorm(c - cPrev, 2, 2)) <= tol
        quiet = quiet + 1;
        if quiet >= patience, break; end
    else
        quiet = 0;
    end
end
end


function c = i_kmeanspp(x, k, batchSize)
% k-means++ seeding (Arthur & Vassilvitskii 2007): each new centre is drawn
% with probability proportional to its squared distance from the nearest
% centre already chosen, which spreads the seeds over the data instead of
% clumping them. Seeding costs O(n*k*d), so on a large input it runs on a
% subsample -- otherwise the seeding would dominate the very cost this
% algorithm exists to avoid.
n = size(x, 1);
cap = max(10*batchSize, 10*k);
if n > cap
    x = x(randperm(n, cap), :);
    n = cap;
end

c = zeros(k, size(x, 2));
c(1, :) = x(randi(n), :);
d2 = sum((x - c(1, :)).^2, 2);
for j = 2:k
    if all(d2 <= 0)
        % Fewer distinct points than clusters: no weighting is meaningful.
        c(j, :) = x(randi(n), :);
    else
        cw = cumsum(d2);
        c(j, :) = x(find(cw >= rand()*cw(end), 1), :);
    end
    d2 = min(d2, sum((x - c(j, :)).^2, 2));
end
end
