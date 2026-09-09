function [coeff, score, latent] = e_kpca(X, k, sgm, isds)
%E_KPCA  Kernel PCA.
%
%   [coeff, score, latent] = PKG.E_KPCA(X, k, sgm, isds)
%
%   SCORE is the embedding of X -- that is the output to plot or to feed to
%   another method. COEFF holds the projection coefficients, for mapping
%   NEW data into the same space, and LATENT the eigenvalues.
%
%   COEFF IS NOT AN EMBEDDING. EIGS returns eigenvectors with unit-norm
%   columns, and SCORE = Kc*coeff = coeff*latent, so coeff(:, j) is
%   score(:, j) divided by its own eigenvalue: taking coeff gives every
%   component equal weight and discards the spectrum. RUN.ML_METAVIZ took
%   output 1 at three call sites and stored it as an embedding, which is
%   the same mistake PKG.E_RANDMDS made by dropping its singular values.
%
% how to use coeff of kpca:
% 	[coeff] = pkg.e_kpca(XTaining);
%   Ktest = constructKernel(XTest)
%   Y = Ktest*coeff;
%   Y is the embedding result of XTest.
%
% (c) 2023 James Cai jcai@tamu.edu

if nargin < 2, k = 2; end
if nargin < 3, sgm = 40; end
if nargin < 4, isds = false; end % if true, then input X is a distance matrix

n = size(X, 1);
if isds
    D = X;
else
    try
        d = dot(X, X, 2);
        D = d + d' - 2 * (X * X');
    catch
        D = pdist2(X, X).^2;
    end
end
K = exp(-D/(2 * sgm^2));
% Centering
H = eye(n) - ones(n, n) / n;
Kc = H * K * H;

%{
https://lvdmaaten.github.io/drtoolbox/
ell = size(X,1);
% Normalize kernel matrix K
mapping.column_sums = sum(K) / ell;                       % column sums
mapping.total_sum   = sum(mapping.column_sums) / ell;     % total sum
J = ones(ell, 1) * mapping.column_sums;                   % column sums (in matrix)
K = K - J - J';
K = K + mapping.total_sum;
%}
[coeff, latent] = eigs(Kc, k);
score = Kc * coeff;
