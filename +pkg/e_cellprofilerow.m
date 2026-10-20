function [difference, level, neighbours] = e_cellprofilerow(model, idx, k)
%E_CELLPROFILEROW  How one cell's genes differ from all cells'.
%
%   [difference, level, neighbours] = pkg.e_cellprofilerow(model, idx, k)
%   takes a MODEL made by PKG.E_CELLPROFILEMODEL and returns, for every
%   gene,
%
%     level       its mean expression over cell IDX and its K-1 nearest
%                 neighbours in principal-component space
%     difference  LEVEL minus the gene's mean over all cells. Expression is
%                 log1p, so a difference of d is roughly a fold change of
%                 exp(d) between the neighbourhood and the average cell.
%
%   NEIGHBOURS holds the K cells averaged, IDX first. K defaults to 15;
%   K = 1 uses the cell alone, which is noisier: most genes are not detected
%   in any one cell.
%
%   See also PKG.E_CELLPROFILEMODEL, GUI.CALLBACK_CELLGENEPROFILE.

arguments
    model (1, 1) struct
    idx (1, 1) double {mustBePositive, mustBeInteger}
    k (1, 1) double {mustBePositive, mustBeInteger} = 15
end

numCells = size(model.Expression, 2);
if idx > numCells
    error("e_cellprofilerow:IndexOutOfRange", ...
        "Cell index %d exceeds the %d cells in the model.", idx, numCells);
end
k = min(k, numCells);

neighbours = knnsearch(model.Searcher, model.Searcher.X(idx, :), K=k).';
% Ties in distance can put another cell at the same point first.
neighbours = [idx; neighbours(neighbours ~= idx)];
neighbours = neighbours(1:k);

level = full(mean(model.Expression(:, neighbours), 2));
difference = level - model.Mean;
end
