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

[ND1 ND2]=size(D);
if K >= ND1
	error('First NN is point itself and hence discarded. Hence K < N, strictly')
end



if (ND1 ~= ND2)
error('Distance matrix is not square!')
end

if(D ~= D')
error('Distance matrix is not symmetric')
end

sorted_dist_mat = zeros(size(D));
sorted_dist_mat=sort(D);

%~ First NN is point itself and discarded. K < N strictly

%disp('ASSUMED THAT DISTANCE MATRIX IS SYMMETRIC')
knnlength = 0;
knnlength = sum(sum(sorted_dist_mat(2:K+1,:)));
return
%~ for i = 1 : ND1
	%~ for j = 1 : K
		%~ knnlength = knnlength + sorted_dist_mat(j+1,i);
	%~ end
%~ end

