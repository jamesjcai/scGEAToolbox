function model = e_cellneighbourmodel(X, g, s, opts)
%E_CELLNEIGHBOURMODEL  What PKG.E_CELLNEIGHBOURROW needs to measure nearness to a cell.
%
%   model = pkg.e_cellneighbourmodel(X, g, s) takes the raw count matrix X
%   (genes x cells), its gene names G and the cell embedding S (cells x d)
%   and returns a struct with fields
%
%     Score          cells x npc principal-component scores (PKG.E_CELLPCS)
%     Embedding      S, as given
%     Transition     cells x cells sparse, row-stochastic one-step transition
%                    matrix of a random walk on the cells' kNN graph
%     NumNeighbours  the k the graph was built with
%
%   model = pkg.e_cellneighbourmodel(..., NumNeighbours=K, NumComponents=NPC)
%   sets the graph's k (default 15) and the number of components (default 30).
%
%   The graph is built in PC space, not on S: t-SNE and UMAP distort the
%   distances between cells, and a walk on them would inherit that. Edge
%   weights are an adaptive Gaussian kernel, exp(-d_ij^2 / (sigma_i*sigma_j))
%   with sigma_i the distance from cell i to its k-th neighbour, symmetrised
%   and row-normalised. Each cell keeps a self-loop of weight 1, so the walk
%   is lazy and does not oscillate between two halves of a bipartite graph.
%
%   See also PKG.E_CELLNEIGHBOURROW, GUI.CALLBACK_CELLNEIGHBOURHOOD.

arguments
    X {mustBeNumeric}
    g
    s {mustBeNumeric}
    opts.NumNeighbours (1, 1) double {mustBePositive, mustBeInteger} = 15
    opts.NumComponents (1, 1) double {mustBePositive, mustBeInteger} = 30
end

numCells = size(X, 2);
if size(s, 1) ~= numCells
    error("e_cellneighbourmodel:SizeMismatch", ...
        "S has %d rows but X has %d cells. Pass one embedding row per cell.", ...
        size(s, 1), numCells);
end
k = min(opts.NumNeighbours, numCells - 1);

score = pkg.e_cellpcs(X, g, opts.NumComponents);

% Column 1 of the search result is each cell itself, at distance 0.
[idx, d] = knnsearch(score, score, K=k + 1);
sigma = d(:, end);
sigma(sigma == 0) = eps;
rows = repmat((1:numCells)', 1, k + 1);
w = exp(-(d.^2) ./ (sigma .* sigma(idx)));
W = sparse(rows(:), idx(:), w(:), numCells, numCells);
W = (W + W.') / 2;
T = W ./ sum(W, 2);

model = struct("Score", score, "Embedding", s, "Transition", T, ...
    "NumNeighbours", k);
end
