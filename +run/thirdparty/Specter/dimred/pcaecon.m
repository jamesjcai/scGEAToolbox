function score = pcaecon(A, k)
% FPCA Fast PCA
% 
%   [U,S,V] = FPCA(A,k,i,usePowerMethod) computes the truncated PCA
%   of the input matrix A upto rank k using i levels of
%   Krylov method as given in [1], p. 3.
% 
%   If usePowerMethod is given as true, then only exponent i is used (i.e.
%   as power method). See [2] p.9, Randomized PCA algorithm for details.
% 
%   [1] Halko, N., Martinsson, P. G., Shkolnisky, Y., & Tygert, M. (2010).
%   An algorithm for the principal component analysis of large data sets.
%   Arxiv preprint arXiv:1007.5510, 0526. Retrieved April 1, 2011, from
%   http://arxiv.org/abs/1007.5510. 
%   
%   [2] Halko, N., Martinsson, P. G., & Tropp, J. A. (2009). Finding
%   structure with randomness: Probabilistic algorithms for constructing
%   approximate matrix decompositions. Arxiv preprint arXiv:0909.4061.
%   Retrieved April 1, 2011, from http://arxiv.org/abs/0909.4061.
% 
%   See also SVD.
% 
%   Copyright 2011 Ismail Ari, http://ismailari.com.

% De-mean
A = bsxfun(@minus,A,mean(A));
%[U, S, V] = fsvd(A, k, i, usePowerMethod);
%score = U(:,1:k)*S(1:k,1:k);

[m,n]=size(A);
[U,S,V] = svdecon(A);
score = U(:,1:k)*S(1:k,1:k);
end
