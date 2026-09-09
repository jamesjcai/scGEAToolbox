function [R] = e_distcorrmtx(X)
%E_DISTCORRMTX Pairwise distance correlation between the rows of X.
%   R = pkg.e_distcorrmtx(X) returns an n-by-n symmetric matrix for the
%   n rows of X, computed the fast way: the doubly-centred distance
%   matrix of each row is formed once, and every pair is then an inner
%   product. It agrees with the pairwise pkg.e_distcorr reference in the
%   comment below to 8e-17.
%
%   THE DIAGONAL IS LEFT AT ZERO, not 1. dCor(x, x) is 1 by definition,
%   so R is an adjacency matrix without self-loops rather than a
%   correlation matrix -- which is what net.distcorrnet wants. Anything
%   that treats R(:) as data has to drop the diagonal first: sc_distcorr
%   regressed dCor on Pearson r without doing so, and the n impossible
%   (r = 1, dCor = 0) points dominated the fit.
%
%   See also NET.DISTCORRNET, PKG.E_DISTCORR, SC_DISTCORR.

[n, m] = size(X);
m2 = m * m;

% Precompute centered distance matrices and variances
dvar = zeros(n, 1);
Avec = cell(n, 1);
for k = 1:n
    xk = X(k, :)';
    a = abs(xk - xk');
    a = a - mean(a) - mean(a, 2) + mean(a(:));
    dvar(k) = a(:)' * a(:) / m2;
    Avec{k} = a(:);
end

% Pairwise distance correlations — cache outer vector only
R = zeros(n);
for k = 1:n - 1
    ak = Avec{k};
    Avec{k} = [];  % free memory as we go
    for l = k + 1:n
        dcov = (ak' * Avec{l}) / m2;
        R(k, l) = sqrt(dcov / sqrt(dvar(k) * dvar(l)));
        R(l, k) = R(k, l);
    end
end
%{
any way to make this code faster:

function [R] = e_distcorrmtx(X)

n = size(X, 1); % number of genes
R = zeros(n);
c = 0;
for k = 1:n - 1
    fprintf('%d......%d\n', k, n)
    for l = k + 1:n
        c = c + 1;
        R(k, l) = pkg.e_distcorr(X(k, :)', X(l, :)');
        R(l, k) = R(k, l);
        % R(k,l)=pkg.e_bcdistcorr(X(k,:)',X(l,:)');
        % if c>20, break; end
    end
    % if c>20, break; end
end
%}