function v = e_cellneighbourrow(model, idx, metric, steps)
%E_CELLNEIGHBOURROW  Every cell's nearness to one focused cell.
%
%   v = pkg.e_cellneighbourrow(model, idx, metric) returns a cells x 1
%   vector relating each cell to cell IDX, from a MODEL made by
%   PKG.E_CELLNEIGHBOURMODEL. METRIC is one of
%
%     "walk"        probability that a random walk of STEPS steps on the
%                   kNN graph, started at cell IDX, ends at each cell. The
%                   values sum to 1. Larger means nearer.
%     "expression"  Euclidean distance in principal-component space.
%                   Smaller means nearer.
%     "embedding"   Euclidean distance in the cell embedding. Smaller means
%                   nearer, as far as t-SNE or UMAP preserve distance.
%
%   v = pkg.e_cellneighbourrow(model, idx, "walk", STEPS) sets the walk
%   length (default 5). Short walks stay in the cell's own neighbourhood;
%   long ones spread along the paths that lead out of it.
%
%   Only one row is computed: the walk is STEPS sparse vector-matrix
%   products, never the dense cells x cells matrix.
%
%   See also PKG.E_CELLNEIGHBOURMODEL, GUI.CALLBACK_CELLNEIGHBOURHOOD.

arguments
    model (1, 1) struct
    idx (1, 1) double {mustBePositive, mustBeInteger}
    metric (1, 1) string {mustBeMember(metric, ["walk", "expression", "embedding"])}
    steps (1, 1) double {mustBePositive, mustBeInteger} = 5
end

numCells = size(model.Transition, 1);
if idx > numCells
    error("e_cellneighbourrow:IndexOutOfRange", ...
        "Cell index %d exceeds the %d cells in the model.", idx, numCells);
end

switch metric
    case "walk"
        p = zeros(1, numCells);
        p(idx) = 1;
        for t = 1:steps
            p = p * model.Transition;
        end
        v = full(p.');
    case "expression"
        v = vecnorm(model.Score - model.Score(idx, :), 2, 2);
    case "embedding"
        v = vecnorm(model.Embedding - model.Embedding(idx, :), 2, 2);
    otherwise
        % The arguments block admits only the three metrics above.
        error("e_cellneighbourrow:UnknownMetric", "Unknown metric ""%s"".", metric);
end
end
