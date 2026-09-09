function [idx, Z] = e_corrmatorder(M)
%E_CORRMATORDER Order a correlation matrix so its blocks sit together.
%   idx = pkg.e_corrmatorder(M) returns a permutation of 1:n that brings
%   correlated variables next to each other, by clustering on 1 - M --
%   the matrix being displayed -- with average linkage and then choosing
%   the leaf order that minimises the distance between neighbours.
%
%   [idx, Z] = pkg.e_corrmatorder(M) also returns the linkage tree.
%
%   GUI.CALLBACK_CELLSCORECORRMATRIX used to reorder its Spearman
%   heatmap with clusterdata(C', maxclust=5) on the raw cell-score
%   matrix. That is Euclidean distance and single linkage on
%   unstandardised scores, so it sorts by each program's score LEVEL,
%   which a rank correlation is invariant to. Two programs with rho = 1
%   on different score scales are adjacent in M and far apart under that
%   metric. Measured on twelve programs in three co-varying blocks whose
%   score levels span 0.01 to 102: the old ordering left 55% of adjacent
%   pairs within a block and interleaved them as 1 2 3 2 2 2 1 1 1 3 3 3;
%   this one reaches 82%, which is every pair except the two block
%   boundaries, as 1 1 1 1 2 2 2 2 3 3 3 3.
%
%   See also PKG.E_CELLSCORECORRMAT, GUI.CALLBACK_CELLSCORECORRMATRIX.

arguments
    M (:, :) double
end

n = size(M, 1);
if n ~= size(M, 2)
    error('pkg:e_corrmatorder:NotSquare', ...
        'M must be square; it is %d-by-%d.', n, size(M, 2));
end

idx = (1:n)';
Z = [];
if n < 3
    % Nothing to order, and LINKAGE needs at least two observations.
    return
end

% 1 - rho, symmetrised against round-off and with an exact zero diagonal,
% which SQUAREFORM requires.
D = 1 - M;
D = (D + D')/2;
D(1:n+1:end) = 0;

if ~all(isfinite(D(:)))
    % A program with no variance gives NaN correlations. Ordering on NaN
    % is meaningless, so keep the original order rather than inventing
    % one; the heatmap still shows the NaNs.
    warning('pkg:e_corrmatorder:NonFiniteEntries', ...
        ['%d of the %d correlations are not finite, so the original ', ...
        'order is kept. A gene program that scores identically in ', ...
        'every cell has no correlation with anything.'], ...
        nnz(~isfinite(D)), numel(D));
    return
end

v = squareform(D, 'tovector');
Z = linkage(v, 'average');
idx = optimalleaforder(Z, v)';
end
