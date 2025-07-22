% [B,D2] = eabounds(logK,D2) Gaussian EAs: bounds
%
% Computes simple bounds (in constant time) of the beta values given the
% distances of each point to its k nearest neighbors, contained in matrix D2:
%   D2(n,i) = squared distance from point X(n,:) to its ith nearest neighbor.
%
% We return D2 because the bounds' formula needs strictly increasing neighbor
% distances, so in case of tied distances we perturb them a little. Tied
% distances can happen with quantized data values, e.g. if the dataset is an
% image that has areas with constant color.
%
% In:
%   logK: scalar, log of the perplexity.
%   D2: N x k matrix of sorted square distances to the k nearest neighbors.
% Out:
%   B: Nx2 matrix of log-beta bounds for each data point.
%   D2: same as D2 with a tiny perturbation of the distance to the first
%      nearest neighbor if K or more nearest neighbors are at exactly the
%      same distance.

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan
function [B, D2] = eabounds(logK, D2)

N = size(D2, 2); % We call "N" the number of nearest neighbors, as in the paper
logN = log(N);
logNK = logN - logK;

delta2 = D2(:, 2) - D2(:, 1);
% Ensure delta2 >= eps, [which means delta2!=0.---AN]
ind = find(delta2 < eps);
i = 3;
flag = 1;
while ~isempty(ind)
    if i > exp(logK) && flag
        % Permute ever so slightly the distance to the first neighbor for points
        % that have K or more closest neighbors at the same distance. Practically
        % this happens only when K is very small.
        D2(ind, 1) = D2(ind, 1) * 0.99;
        flag = 0;
    end
    delta2(ind) = D2(ind, i) - D2(ind, 1);
    ind = find(delta2 < eps);
    i = i + 1;
end

deltaN = D2(:, N) - D2(:, 1);

% Compute p1(N,logK)
if logK > log(sqrt(2*N))
    p1 = 3 / 4;
else
    p1 = 1 / 4;
    for i = 1:100, e = -p1 * log(p1/N) - logK;
        g = -log(p1/N) - 1;
        p1 = p1 - e / g;
    end
    p1 = 1 - p1 / 2;
end

bU1 = (2 * log(p1*(N - 1)/(1 - p1))) ./ delta2;
bL1 = (2 * logNK / (1 - 1 / N)) ./ deltaN;
bL2 = (2 * sqrt(logNK)) ./ sqrt(D2(:, N).^2-D2(:, 1).^2);
B = log([max(bL1, bL2), bU1]);