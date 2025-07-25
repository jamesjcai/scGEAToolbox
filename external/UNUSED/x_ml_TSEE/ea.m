% [W,s] = ea(X,K[,k]) Gaussian entropic affinities
%
% This computes the matrix W of Gaussian entropic affinities (EAs) and the
% corresponding bandwidth values s for a dataset X and a desired perplexity
% K. By default, W is a full matrix. Use the optional argument k (with k > K)
% to restrict the EAs to the k nearest neighbors of each point, in which case
% W will be sparse. Note that, if using k neighbors, the maximum perplexity
% achievable is K=k, in the limit where s=Inf (which results in affinities
% equal to 1/k for each neighbor), so k < K makes no sense. In practice, K
% should not approach k too much.
%
% W is normalized so each row sums 1, thus it is a stochastic matrix. It is
% the random-walk matrix commonly used in machine learning, but where each
% row (data point) has its own bandwidth in the Gaussian kernel.
%
% When using sparse affinities, the runtime is mostly due to computing the
% nearest neighbors in nnsqdist.m, which is O(D.N?, rather than to computing
% the EAs, which is O(N). If you have precomputed the nearest neighbors, you
% should pass them as input argument to ea.m. See also nnsqdist.m about using
% selection rather than sorting.
%
% In:
%   X: NxD matrix of N row D-dimensional vectors.
%   K: scalar in (0,N), the perplexity.
%   k: either a scalar k > K, the number of neighbors; or a cell array {D2,nn}
%      containing two N x k matrices, where D2(n,i) is the square distance to
%      the ith nearest neighbor and nn(n,i) is the index of the ith nearest
%      neighbor. Default: N-1.
% Out:
%   W: NxN matrix of EAs (random-walk matrix).
%   s: Nx1 vector of bandwidth values.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [W, s] = ea(X, K, k)

N = size(X, 1);
% ---------- Argument defaults ----------
if ~exist('k', 'var') || isempty(k), k = N - 1; end
if iscell(k)
    D2 = k{1};
    nn = k{2};
    k = size(D2, 2);
else
    [D2, nn] = nnsqdist(X, k); % Square distances to k nn
    %[D2,nn] = nnsqdist(X,k,'mink');	% Usually faster, see nnsqdist.m
end
% ---------- End of "argument defaults" ----------

b = zeros(N, 1);
Wp = zeros(N, k);
logK = log(K);
[B, D2] = eabounds(logK, D2); % Log-beta bounds
[~, p] = sort(D2(:, ceil(K))); % Point order: distance to Kth nn
j = p(1);
b0 = mean(B(j, :));
p = [p; 0]; % Initialization
for i = 1:N % Compute log-beta & EAs for each point
    [b(j), Wp(j, :)] = eabeta(D2(j, :), b0, logK, B(j, :));
    b0 = b(j);
    j = p(i+1); % Next point
end
W = sparse(repmat((1:N)', 1, k), nn, Wp, N, N);
if k >= N - 1, W = full(W);
end
if nargout == 2, s = 1 ./ sqrt(2*exp(b)); end % Bandwidths from log-beta values
