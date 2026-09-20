function [gid, ngroups] = i_partitioncells(Y, opts)
%I_PARTITIONCELLS Cut cells into compact groups of roughly equal size.
%
%   gid = pkg.i_partitioncells(Y)
%   gid = pkg.i_partitioncells(Y, TargetSize=50)
%   [gid, ngroups] = pkg.i_partitioncells(___)
%
%   Splits the rows of Y recursively at the median of their own first
%   principal direction, stopping when a piece is small enough (principal
%   direction divisive partitioning, Boley 1998). Every split is exactly
%   balanced, so the leaves come out within a factor of two of each other
%   in size - which is the point. The groups are pooling neighbourhoods,
%   not cell types: something downstream is going to spend a fixed budget
%   of work per group, and that budget only buys the same precision
%   everywhere if the groups are the same size.
%
%   This is not a substitute for clustering and should not be reported as
%   one. It cuts through whatever is there, including through the middle
%   of a homogeneous population, and a group is only useful on the
%   assumption that a few dozen mutual neighbours in Y share a cell type.
%   Whether they do is a property of Y, so pass a space in which that is
%   plausible - PCA coordinates of log-normalised counts, not a 2-D UMAP,
%   whose distortions put unrelated cells next to each other.
%
%   WHY NOT K-MEANS WITH K = N/TARGETSIZE. k-means gives groups whose
%   sizes vary by an order of magnitude, so a fixed per-group budget
%   oversamples the small ones and starves the large ones; and k-means++
%   seeding is O(N*K), which at K = N/50 is quadratic in the cell count.
%   This is O(N*D^2*log(N/TargetSize)) and needs no seeding, so it is also
%   deterministic: the same Y always gives the same partition.
%
%   INPUTS:
%     Y          - NumCells-by-D coordinates, one cell per row. D is meant
%                  to be small (tens); the per-node work is quadratic in
%                  it, because the splitting direction comes from an
%                  eigendecomposition of the node's D-by-D scatter matrix.
%     TargetSize - (50) wanted group size. Leaves land between
%                  0.75*TargetSize and 1.5*TargetSize, averaging near it;
%                  exact equality is not achievable by bisection and is
%                  not worth chasing.
%
%   OUTPUTS:
%     gid     - NumCells-by-1 group index in 1:ngroups, numbered by the
%               first cell each group contains
%     ngroups - number of groups
%
% See also PKG.I_POOLEDANNOTATE, PKG.I_MAJORITYVOTE, SC_CLUSTER_S.

arguments
    Y (:,:) {mustBeNumeric, mustBeReal, mustBeNonempty}
    opts.TargetSize (1,1) double {mustBeInteger, mustBePositive} = 50
end

Y = double(Y);
if any(~isfinite(Y), 'all')
    error('pkg:i_partitioncells:nonFinite', ...
        ['Y contains Inf or NaN, which have no median and no principal ', ...
        'direction. Drop or impute those cells before partitioning.']);
end

n = size(Y, 1);
gid = ones(n, 1);
ngroups = 1;

% A node is split while it is larger than this, so a leaf is at most
% MAXLEAF and - because every split is exactly balanced - at least half
% of it. Put the cut at 1.5*TargetSize and leaves straddle TargetSize
% instead of sitting above it.
maxleaf = ceil(1.5*opts.TargetSize);
if n <= maxleaf, return; end

% Both stacks are bounded: a balanced bisection down to MAXLEAF/2 leaves
% fewer than 2*n/TargetSize pieces, and the frontier is smaller still.
capacity = 2*ceil(n/opts.TargetSize) + 2;
stack = cell(capacity, 1);
leaves = cell(capacity, 1);
nstack = 1;
nleaves = 0;
stack{1} = (1:n).';

while nstack > 0
    rows = stack{nstack};
    nstack = nstack - 1;
    if numel(rows) <= maxleaf
        nleaves = nleaves + 1;
        leaves{nleaves} = rows;
        continue;
    end
    [left, right] = i_bisect(Y(rows, :), rows);
    stack{nstack+1} = left;
    stack{nstack+2} = right;
    nstack = nstack + 2;
end

% Number the groups by the first cell in each, so GID is a function of the
% partition alone and not of the order the stack happened to unwind in.
leaves = leaves(1:nleaves);
firstcell = cellfun(@min, leaves);
[~, order] = sort(firstcell);
for k = 1:nleaves
    gid(leaves{order(k)}) = k;
end
ngroups = nleaves;
end


function [left, right] = i_bisect(Yn, rows)
% Split ROWS in half at the median of their first principal direction.
%
% The halves are taken from the sorted projection rather than by comparing
% against MEDIAN(p): with many cells tied at the median - duplicate rows,
% or a coordinate that is constant over the node - a `p <= median(p)` test
% can put every cell on one side and loop forever. SORT is stable, so ties
% fall in cell order and the split stays deterministic.

Yn = Yn - mean(Yn, 1);
scattermat = Yn.'*Yn;
% Exact symmetry, so EIG takes the symmetric path and returns real
% eigenvectors in a fixed order. Not named SCATTER: that shadows the
% plotting function for the rest of the file.
scattermat = (scattermat + scattermat.')/2;
[V, d] = eig(scattermat, 'vector');
[~, j] = max(d);
projection = Yn*V(:, j);

[~, order] = sort(projection);
half = floor(numel(rows)/2);
left = rows(order(1:half));
right = rows(order(half+1:end));
end
