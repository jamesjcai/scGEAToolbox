function A = pearsonnet(X, corr_thr)
% Construct co-expression network from thresholded Pearson correlation
%
% A = net.pearsonnet(X)
% A = net.pearsonnet(X, corr_thr)
%
% X        - genes x cells expression matrix (log-normalised recommended)
% corr_thr - absolute correlation threshold to keep edges (default: 0.3)
% A        - genes x genes sparse adjacency matrix
%
% Uses only base MATLAB (no Statistics Toolbox required).

arguments
    X double
    corr_thr(1,1) double {mustBeInRange(corr_thr, 0, 1)} = 0.3
end

X = double(full(X));

% Z-score each gene across cells
mu  = mean(X, 2);
sig = std(X, 0, 2);
sig(sig < eps) = 1;
Xz  = (X - mu) ./ sig;

nc = size(X, 2);
C  = (Xz * Xz') ./ (nc - 1);   % Pearson correlation matrix
C  = max(min(C, 1), -1);

% Zero out weak and self correlations, then sparsify
C(abs(C) < corr_thr) = 0;
C(1:size(C,1)+1:end) = 0;
A = sparse(C);
end
