function [idx, C, D, DX] = diffusion_kmeans(X, k, phi0, Niter, epsilon)

% DIFFUSION_KMEANS diffusion K-means clustering
% function [idx, C, D, DX] = diffusion_kmeans(X, k, phi0, Niter, epsilon)
%
% Input:
%   X(n,d)     --- diffusion coords (right eigenvecs rescaled with
%                  eigenvals, skip first trivial one); each row is
%                  a data pt  
%   k          --- number of clusters
%   phi0       --- first left eigenvector of P (prop. to stationary distr)
%   Niter      --- number of iterations to repeat the clustering with new IC
%   epsilon    --- relative distortion (stopping criteria); default: 0.001
%
% Output:
%   idx(n,1)   --- labeling of data points; labels between 1 and k
%   C(k,d)     --- geometric centroids
%   D          --- distortion (sum of squared distances)
%   DX(n,k)    --- squared distance from each point to every centroid
%
% Calls:       --- distortionMinimization (in each iteration)
%
% Ann B. Lee, May 2008. Last changed, JWR: 3/23/2009

[N,d]=size(X);
aD=Inf; % maximum distortion
if(nargin<5)
    epsilon=1e-3;
end

%--------------------------------------
% K-MEANS (repeat Niter times)
%--------------------------------------
for iter=1:Niter
    tmp_ind = ceil(rand(k,1)*N);  % k random points in X
    c_0 = X(tmp_ind,:);   % k-by-d matrix of initial centroids
    [idx,c,cindex,D,DX]=distortionMinimization(X,phi0,k,c_0,0,epsilon);
    if(D<aD) % keep best result
        aD=D; aDX=DX;
        a_idx=idx;
        ac=c;
    end
end
D=aD; DX=aDX; idx=a_idx; C=ac;



