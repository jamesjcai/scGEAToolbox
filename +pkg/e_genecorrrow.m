function rho = e_genecorrrow(model, idx, method)
%E_GENECORRROW  Correlation of one gene with every gene, across cells.
%
%   rho = pkg.e_genecorrrow(model, idx) returns a genes x 1 vector of the
%   Pearson correlation, across cells, of log1p normalised expression
%   between gene IDX and each gene, from a MODEL made by
%   PKG.E_GENECORRMODEL. A gene that is constant across cells has no
%   correlation, and its entry is NaN.
%
%   rho = pkg.e_genecorrrow(model, idx, "spearman") uses the Spearman
%   correlation instead. The model must have been made with Spearman=true.
%
%   Only one row is computed, with one sparse matrix-vector product; the
%   genes x genes matrix never exists.
%
%   See also PKG.E_GENECORRMODEL, GUI.CALLBACK_GENECORRPARTNERS.

arguments
    model (1, 1) struct
    idx (1, 1) double {mustBePositive, mustBeInteger}
    method (1, 1) string {mustBeMember(method, ["pearson", "spearman"])} = "pearson"
end

numGenes = numel(model.Genes);
if idx > numGenes
    error("e_genecorrrow:IndexOutOfRange", ...
        "Gene index %d exceeds the %d genes in the model.", idx, numGenes);
end

if method == "spearman"
    if isempty(model.Spearman)
        error("e_genecorrrow:NoRanks", ['The model holds no ranks. ' ...
            'Make it with pkg.e_genecorrmodel(X, g, Spearman=true).']);
    end
    m = model.Spearman;
else
    m = model.Pearson;
end
n = size(m.Matrix, 2);
y = m.Matrix(idx, :).';
covariance = (full(m.Matrix*y) - n*m.Mean*m.Mean(idx))/(n - 1);
rho = covariance./(m.Std*m.Std(idx));
% Rounding can put a gene's correlation with itself a hair past 1. Clamp
% before marking the constant genes: MAX and MIN skip NaN, so clamping
% after would turn every NaN into -1.
rho = min(max(rho, -1), 1);
rho(m.Std == 0 | m.Std(idx) == 0 | ~isfinite(covariance)) = NaN;
end
