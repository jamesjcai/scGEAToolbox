function idx = sc_geosketch(X, N, options)
%SC_GEOSKETCH Geometric sketching subsample of cells.
%   idx = SC_GEOSKETCH(X, N) selects N cells whose low-dimensional embedding
%   evenly covers the geometry of the data, over-representing rare cell
%   populations relative to uniform random sampling.
%
%   This is a native MATLAB reimplementation of geosketch
%   (Hie et al., Cell Systems 2019, PMID:31176620). The data geometry is
%   covered with a grid of equal-sized boxes whose resolution is chosen so
%   that at least N boxes are occupied; one cell is then drawn uniformly at
%   random from each of N randomly chosen occupied boxes.
%
%   Inputs:
%     X : numCells-by-numDims embedding (e.g., PCA scores). Rows are cells,
%         columns are dimensions.
%     N : number of cells to sample.
%
%   Name-value arguments:
%     Seed : RNG seed for reproducibility (default [], no reseed).
%
%   Output:
%     idx : N-by-1 indices of the sampled cells (into the rows of X).
%
%   Example:
%     [~, sc] = pca(log1p(sc_norm(sce.X))', "NumComponents", 100);
%     idx = sc_geosketch(sc, 2000);
%     sce = sce.selectcells(idx);

arguments
    X (:, :) {mustBeNumeric, mustBeReal}
    N (1, 1) {mustBePositive, mustBeInteger}
    options.Seed = []
end

numCells = size(X, 1);
if N >= numCells
    idx = (1:numCells)';
    return;
end

if ~isempty(options.Seed)
    rng(options.Seed);
end

X = double(X);

% Scale each dimension to [0, 1] so the covering boxes are isotropic.
xmin = min(X, [], 1);
xmax = max(X, [], 1);
xrange = xmax - xmin;
xrange(xrange == 0) = 1;            % guard against constant dimensions
Xs = (X - xmin) ./ xrange;

% Find the coarsest grid whose covering has at least N occupied boxes.
[boxid, numBoxes] = in_covering(Xs, N);

% Sample distinct occupied boxes uniformly, then one random cell per box.
boxSample = randperm(numBoxes, min(N, numBoxes));
idx = zeros(numel(boxSample), 1);
for k = 1:numel(boxSample)
    members = find(boxid == boxSample(k));
    idx(k) = members(randi(numel(members)));
end

% If the covering yielded fewer boxes than N (heavily duplicated
% coordinates), top up with random unpicked cells to return exactly N.
if numel(idx) < N
    remaining = setdiff((1:numCells)', idx);
    need = N - numel(idx);
    extra = remaining(randperm(numel(remaining), need));
    idx = [idx; extra];
end
idx = idx(:);
end


function [boxid, numBoxes] = in_covering(Xs, target)
% Assign each cell an occupied-box label for the grid whose number of
% occupied boxes is as close as possible to (and at least) target. The box
% count must land near target: a far finer grid drives every cell into its
% own box, which degrades the sketch back to uniform sampling. A doubling
% search brackets the resolution, then a bisection on a continuous number of
% bins per dimension homes in on the smallest resolution reaching target.
% Occupied boxes are the unique rows of the bin-index matrix, so the box
% count is bounded by the number of cells.
numCells = size(Xs, 1);

% Doubling search: grow bins/dimension until the target box count is reached.
lo = 1;
hi = 1;
[boxid, numBoxes] = in_occupancy(Xs, hi);
while numBoxes < target && hi < numCells
    lo = hi;
    hi = hi*2;
    [boxid, numBoxes] = in_occupancy(Xs, hi);
end
if numBoxes < target
    return;             % finest useful grid still below target; use as-is
end

% Bisection for the smallest (continuous) bins/dimension reaching target.
for iter = 1:40
    mid = (lo + hi)/2;
    [idMid, nbMid] = in_occupancy(Xs, mid);
    if nbMid >= target
        hi = mid;
        boxid = idMid;
        numBoxes = nbMid;
    else
        lo = mid;
    end
    if hi - lo < 1e-6, break; end
end
end


function [id, nb] = in_occupancy(Xs, nbins)
% Bin the scaled coordinates and label each distinct occupied box.
B = floor(Xs*nbins);
B = min(B, nbins - 1);             % keep the maximum value inside the last bin
[~, ~, id] = unique(B, "rows");
nb = max(id);
end
