function [E,v,score,K] = km_kpca(X,m,ktype,kpar)
% KM_KPCA performs kernel principal component analysis (KPCA) on a data set
% X.
% Input:	- X: data matrix in column format (each data point is a row)
%			- m: the number of principal components to return. If m is 
%			smaller than 1, it is interpreted as the fraction of the signal
%			energy that is to be contained within the returned principal
%			components.
%			- ktype: string representing kernel type.
%			- kpar: vector containing the kernel parameters.
% Output:	- E: matrix containing the principal components.
%			- v: array containing the eigenvalues.
% USAGE: [E,v] = km_kpca(X,m,ktype,kpar)
%
% Author: Steven Van Vaerenbergh (steven.vanvaerenbergh at unican.es) 2010.
%
% This file is part of the Kernel Methods Toolbox for MATLAB.
% https://github.com/steven2358/kmbox

if nargin<4, kpar=1.0; end
if nargin<3, ktype='gauss'; end

n = size(X,1);
K = pkg.km_kernel(X,X,ktype,kpar);


H = eye(n) - ones(n,n)/n;
K_centered = H * K * H;

[E,V] = eigs(K_centered,m);
v=diag(V);
%[E,V] = eig(K);
%v = diag(V); % eigenvalues
%[v,ind] = sort(v,'descend');
%v = v(1:m);
%E = E(:,ind(1:m)); % principal components

score=(E'*K_centered)';

% for k=1:m
%     E(:,k) = E(:,k)/sqrt(n*v(k)); % normalization
% end