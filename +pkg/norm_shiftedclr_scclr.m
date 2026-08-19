function [Xclr, alpha] = norm_shiftedclr_scclr(X, alpha)
% norm_shiftedclr_scclr  Canonical scclr shifted-CLR (PFlog) normalization.
%   Xclr = norm_shiftedclr_scclr(X) estimates the negative-binomial
%   overdispersion alpha from the data and returns the shifted CLR
%   center(log1p(4*alpha*X)).
%   Xclr = norm_shiftedclr_scclr(X, alpha) uses a caller-supplied alpha
%   instead of estimating it.
%   [Xclr, alpha] = norm_shiftedclr_scclr(...) also returns the alpha used.
%
%   X is genes-by-cells (sparse or dense) raw counts. Xclr is a dense
%   genes-by-cells matrix; per-cell mean subtraction destroys sparsity.
%
%   This matches the current cleartools/scclr definition (target="auto"),
%   where the shifted CLR is center(log(x + 1/(4*alpha))), computed as the
%   equivalent sparsity-preserving center(log1p(4*alpha*x)). Unlike
%   pkg.norm_shiftedclr, which applies a per-cell proportional-fitting scale
%   mean(depth)/depth, this uses a single global scale 4*alpha (the cell
%   depth cancels).
%
%   Reference: cleartools/scclr and cleartools/runorm.
%   https://www.biorxiv.org/content/10.1101/2022.05.06.490859v3

arguments
    X {mustBeNumeric}
    alpha {mustBeScalarOrEmpty} = []
end

% sparse stays sparse (already double); integer/single densify to double
X = double(X);

if isempty(alpha)
    alpha = i_estimate_overdispersion(X);
end

if ~(alpha > 0)
    error("norm_shiftedclr_scclr:NonPositiveOverdispersion", ...
        "Overdispersion alpha = %g is non-positive, so the shift " + ...
        "1/(4*alpha) is undefined. Supply a positive alpha, e.g. " + ...
        "norm_shiftedclr_scclr(X, 0.5).", alpha);
end

scale = 4 * alpha;

% log1p(4*alpha*x) on nonzeros; zeros stay zero because log1p(0) = 0
Xs = X;
Xs(Xs~=0) = log1p(scale * Xs(Xs~=0));

% CLR centering: subtract per-cell mean over all genes, including zeros
cell_mean = mean(Xs, 1);

Xclr = full(Xs) - cell_mean;

end

function alpha = i_estimate_overdispersion(X)
% Estimate NB overdispersion from Var_g = mu_g + alpha*mu_g^2 across genes.
% The model is linear in alpha, giving the closed-form OLS solution
% alpha = sum_g (Var_g - mu_g)*mu_g^2 / sum_g mu_g^4. X is genes-by-cells,
% so per-gene moments are taken over cells (dim 2) with divisor n (the
% population variance).
n = size(X, 2);
mu = full(sum(X, 2)) / n;
ex2 = full(sum(X.^2, 2)) / n;
variance = ex2 - mu.^2;
mu2 = mu.^2;

num = sum((variance - mu) .* mu2);
den = sum(mu2.^2);

if den == 0
    error("norm_shiftedclr_scclr:AllZero", ...
        "All genes are identically zero, so overdispersion is undefined.");
end

alpha = num / den;

end
