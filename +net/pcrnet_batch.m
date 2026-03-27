function [A, coeffv] = pcrnet_batch(X, ncom)
% Construct GRN using principal component regression (batch pagesvd)
%
% A = net.pcrnet_batch(X)
% [A, coeffv] = net.pcrnet_batch(X, ncom)
%
% X    - genes x cells expression matrix (LogNormalized recommended)
% ncom - number of principal components (default: 3)
%
% Stacks all leave-one-out matrices into a 3D tensor and uses pagesvd
% to batch-decompose all pages at once, then loops only for the
% lightweight regression step.
%
% ref: https://rdrr.io/cran/dna/man/PCnet.html
%      https://github.com/cran/dna/blob/master/src/rpcnet.c

if nargin < 2 || isempty(ncom), ncom = 3; end

X = X.';
X = zscore(X);
[m, n] = size(X);
A = 1 - eye(n);

% Pre-allocate leave-one-out tensor
D = zeros(m, n - 1, n);
idx_all = 1:n;
for k = 1:n
    D(:, :, k) = X(:, [idx_all(1:k-1), idx_all(k+1:end)]);
end
disp('Step 1.')

[~, ~, coeffv] = pagesvd(D, "econ");
disp('Step 2.')

for k = 1:n
    y = X(:, k);
    Xi = D(:, :, k);
    coeff = coeffv(:, 1:ncom, k);
    score = Xi * coeff;
    nrm2 = sum(score .* score);
    Beta = sum(y .* (score ./ nrm2));
    A(k, A(k, :) == 1) = coeff * Beta';
end
disp('Step 3.')
end
