function [c, info] = sc_louvain(s, k, opts)
% sc_louvain - cluster cells by Louvain community detection on a kNN graph
%
%   C = SC_LOUVAIN(S, K) builds a shared-nearest-neighbour graph from the
%   cell embedding S, runs Louvain on it, and returns a cluster index per
%   cell. The resolution is tuned by bisection so that the number of
%   communities comes out at K.
%
%   C = SC_LOUVAIN(S) uses the default resolution and lets the graph
%   decide how many communities there are.
%
%   C = SC_LOUVAIN(..., Resolution=G) fixes the resolution at G and
%   ignores K. Larger values give more, smaller clusters.
%
%   C = SC_LOUVAIN(..., NumNeighbors=N) sets the neighbourhood size used
%   to build the graph (default 20), and Prune=P drops graph edges whose
%   Jaccard weight is below P (default 1/15).
%
%   [C, INFO] = SC_LOUVAIN(...) returns the resolution used, the number of
%   communities found, the modularity, and the graph.
%
%   Clusters are numbered by size, so cluster 1 is the largest.
%
%   This is the graph-based route the single-cell field settles on: the
%   graph is built the way Seurat's FindNeighbors does it (kNN, then
%   Jaccard-weighted shared neighbours, then pruning), and clustered the
%   way FindClusters does. Clustering the graph rather than the 2-D or
%   3-D embedding also avoids inheriting the distortions t-SNE and UMAP
%   introduce when they flatten the data.
%
%   Example:
%       s = sc_umap(X, 2);
%       c = sc_louvain(s, 8);
%
% see also: sc_cluster_s, pkg.e_louvain, sc_knngraph

arguments
    s {mustBeNumeric}
    k double = []
    opts.NumNeighbors (1, 1) double {mustBePositive} = 20
    opts.Resolution double = []
    opts.Prune (1, 1) double = 1 / 15
    opts.MaxTuningSteps (1, 1) double = 25
end

n = size(s, 1);
if n < 3
    error('sc_louvain:tooFewCells', ...
        'S has %d rows. At least 3 cells are needed to build a graph.', n);
end

% ---- shared-nearest-neighbour graph, as Seurat builds it ---------------
knnK = min(opts.NumNeighbors, n-1) + 1;   % column 1 is the cell itself
nnIdx = knnsearch(s, s, K=knnK);
rowIdx = repmat((1:n)', 1, knnK);
A = sparse(rowIdx, nnIdx, 1, n, n);
A = double(A > 0);

shared = A * A';
[gi, gj, gv] = find(shared);
w = gv ./ (2 * knnK - gv);            % Jaccard index of the two neighbourhoods
keep = gi ~= gj & w >= opts.Prune;
W = sparse(gi(keep), gj(keep), w(keep), n, n);
W = max(W, W');

% ---- resolution ---------------------------------------------------------
if ~isempty(opts.Resolution)
    gammaUsed = opts.Resolution;
    [c, Q] = pkg.e_louvain(W, gammaUsed);
elseif isempty(k) || k <= 1
    gammaUsed = 0.8;
    [c, Q] = pkg.e_louvain(W, gammaUsed);
else
    [c, gammaUsed, Q] = i_tuneresolution(W, k, opts.MaxTuningSteps);
    if max(c) ~= k
        warning('sc_louvain:clusterCountMissed', ...
            ['Asked for %d clusters, found %d at resolution %.4g. ', ...
            'Louvain sets the number of communities by the structure of ', ...
            'the graph, so an exact count is not always reachable; pass ', ...
            'Resolution= to set it directly.'], k, max(c), gammaUsed);
    end
end

% ---- number the clusters by size ---------------------------------------
counts = accumarray(c, 1);
[~, order] = sort(counts, 'descend');
relabel = zeros(numel(counts), 1);
relabel(order) = 1:numel(counts);
c = relabel(c);

if nargout > 1
    info = struct('Resolution', gammaUsed, 'NumClusters', max(c), ...
        'Modularity', Q, 'Graph', W, 'NumNeighbors', opts.NumNeighbors);
end
end

function [cBest, gBest, qBest] = i_tuneresolution(W, k, maxSteps)
% Bisect the resolution towards K communities. The count rises with the
% resolution but in steps, so an exact K is not always on the curve; keep
% the closest partition seen.
cBest = [];
gBest = NaN;
qBest = NaN;
nBest = inf;

    function keepIfBetter(c, g, q)
        nc = max(c);
        if abs(nc-k) < abs(nBest-k) || (abs(nc - k) == abs(nBest - k) && g < gBest)
            cBest = c;
            gBest = g;
            qBest = q;
            nBest = nc;
        end
    end

gLo = 0.05;
gHi = 2;
[cc, qq] = pkg.e_louvain(W, gLo);
keepIfBetter(cc, gLo, qq);
nLo = max(cc);
step = 0;
while nLo > k && step < 8
    gLo = gLo / 4;
    [cc, qq] = pkg.e_louvain(W, gLo);
    keepIfBetter(cc, gLo, qq);
    % A graph with several connected components has a floor on the number
    % of communities that no resolution gets under. Stop when lowering it
    % stops helping, rather than running the bracket down to 1e-7.
    if max(cc) >= nLo, break; end
    nLo = max(cc);
    step = step + 1;
end

[cc, qq] = pkg.e_louvain(W, gHi);
keepIfBetter(cc, gHi, qq);
nHi = max(cc);
step = 0;
while nHi < k && step < 12
    gHi = gHi * 2;
    [cc, qq] = pkg.e_louvain(W, gHi);
    keepIfBetter(cc, gHi, qq);
    if max(cc) <= nHi, break; end
    nHi = max(cc);
    step = step + 1;
end

for it = 1:maxSteps
    if nBest == k || gHi - gLo < 1e-4, break; end
    gMid = (gLo + gHi) / 2;
    [cc, qq] = pkg.e_louvain(W, gMid);
    keepIfBetter(cc, gMid, qq);
    nMid = max(cc);
    if nMid < k
        gLo = gMid;
    elseif nMid > k
        gHi = gMid;
    else
        break
    end
end
end
