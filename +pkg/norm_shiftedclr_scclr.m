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
%   mean(depth)/depth, this uses a single global scale 4*alpha and no size
%   factor at all.
%
%   NOTE ON DEPTH. The header used to assert that depth cancels out here.
%   It does not, quite. The cancellation log(d_c*p_g) = log(d_c) + log(p_g) needs the log
%   applied to all G entries, and two things break it here: LOG1P is only
%   additive for arguments well above 1, and the zeros are deliberately left
%   at 0 rather than shifted, so they carry no -log(d_c) term while the
%   per-cell mean below averages over them anyway. Measured on 400 cells
%   sharing one composition with depths from 21 to 21737, each gene's
%   normalised value still correlated with log depth at a median Spearman
%   rho of +0.100, against -0.015 for PKG.NORM_SHIFTEDCLR, which does apply
%   the proportional-fitting scale. Use that one if depth correction
%   matters; the transform here is left as the scclr definition it names.
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
% Estimate the NB overdispersion alpha, with the cell size factors in the
% model.
%
%   For x_gc ~ NB(mean = s_c*q_g, overdispersion alpha),
%
%       E[x_gc]   = s_c*q_g
%       Var[x_gc] = s_c*q_g + alpha*(s_c*q_g)^2
%
%   Estimating q_g as m_g = sum_c x_gc / sum_c s_c and writing
%   SSE_g = sum_c (x_gc - s_c*m_g)^2 gives
%
%       E[SSE_g] = m_g*S1 + alpha*m_g^2*S2,   S1 = sum s_c, S2 = sum s_c^2
%
%   and the same mu^2-weighted OLS across genes as before then gives
%
%       alpha = sum_g (SSE_g - m_g*S1)*m_g^2 / (S2 * sum_g m_g^4).
%
%   WHY THE SIZE FACTORS ARE NEEDED. This used to take per-gene moments over
%   cells with no size factors at all, treating raw counts as if every cell
%   had the same depth. For a gene at composition p_g in a cell of depth
%   d_c, the variance over cells is mu_g + alpha*mu_g^2 + p_g^2*Var(d);
%   dividing that third term by mu_g^2 leaves Var(d)/E[d]^2, the squared
%   coefficient of variation of the depth, added to the slope for every
%   gene. So the estimate measured depth spread rather than overdispersion.
%
%   Measured on 300 genes and 400 cells drawn from an NB with a known alpha:
%
%       CV depth   true alpha   old estimate   this estimate
%          0.00       0.10          0.0949         0.0914
%          0.32       0.10          0.2132         0.0949
%          0.88       0.10          0.9681         0.0959
%          1.94       0.10          4.3702         0.0870
%          1.94       0.30          5.5367         0.2153
%          1.94       1.00          8.3928         0.7345
%
%   Real scRNA-seq has a depth CV around 0.3 to 1.0, so the old estimate ran
%   2 to 10 times high. That matters here beyond the reported value: alpha
%   sets the shift 1/(4*alpha), so an alpha ten times too large made the
%   shift ten times too small and the transform far more log-like than the
%   method specifies. When every cell has the same depth all s_c are 1 and
%   the expression above reduces to the previous one exactly, which the
%   first row of the table shows.

d = full(sum(X, 1));
if ~any(d > 0)
    error("norm_shiftedclr_scclr:AllZero", ...
        "All genes are identically zero, so overdispersion is undefined.");
end

% Size factors normalised to mean 1 over the cells that have counts, so
% that m_g stays on the scale of a count. An empty cell gets s_c = 0 and
% contributes nothing to S1, S2 or SSE, which is right.
s = d / mean(d(d > 0));
S1 = sum(s);
S2 = sum(s.^2);

m = full(sum(X, 2)) / S1;

% SSE_g = sum_c (x_gc - s_c*m_g)^2, expanded so the residuals are never
% formed and X can stay sparse.
sse = full(sum(X.^2, 2)) - 2*m.*full(X*s(:)) + (m.^2)*S2;

m2 = m.^2;
num = sum((sse - m*S1) .* m2);
den = S2 * sum(m2.^2);

if den == 0
    error("norm_shiftedclr_scclr:AllZero", ...
        "All genes are identically zero, so overdispersion is undefined.");
end

alpha = num / den;

end
