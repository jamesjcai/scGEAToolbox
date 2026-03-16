function knnlength = knn_distmat(D, K)
%KNN K-Nearest Neighbor classifier using an arbitrary distance matrix
%
%  knnlength = knn(D, K)
%
%  Input and output arguments:
%   D     (matrix) of size NxN: This is a precalculated dissimilarity (distance matrix).
% NOTE:  D MUST BE  SYMMETRIC!
%   K    (scalar) the K in K-NN classifier
%   knnlength:   Length of kNN graph
%
% Copyright (c) Huzefa Neemuchwala,( hneemuch@umich.edu ). All rights reserved.

[ND1, ND2] = size(D);
if K >= ND1
    error('First NN is point itself and hence discarded. Hence K < N, strictly')
end

if ND1 ~= ND2
    error('Distance matrix is not square!')
end

if ~issymmetric(D) % ~= D.'
    error('Distance matrix is not symmetric')
end

% Exclude self-distances explicitly so duplicate zeros do not bias K-NN selection.
Doff = reshape(D(~eye(ND1)), ND1 - 1, ND1);
sorted_dist_mat = sort(Doff, 1);

knnlength = sum(sum(sorted_dist_mat(1:K, :)));