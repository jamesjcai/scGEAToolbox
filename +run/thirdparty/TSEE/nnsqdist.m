% [D2,nn] = nnsqdist(X,k[,method]) Nearest-neighbor squared distances
%
% This computes the k nearest neighbors of each data point in X and its
% squared Euclidean distances:
% - D2(n,i) = squared distance from point X(n,:) to its ith nearest neighbor.
% - nn(n,i) = index of the ith nearest neighbor of X(n,:).
%
% Finding the nearest neighbors is done by default by sorting N distances for
% each point. If k<<N, a faster way is to use selection in O(N+k.logk) rather
% than sorting in O(N.logN). To use selection, install this package and set
% method='mink':
%   Min/Max selection
%   http://www.mathworks.com/matlabcentral/fileexchange/23576-minmax-selection
% Its mink() function selects the k smallest elements in the array with a
% partial quicksort in O(N+k.logk).
% Note that computing the distances between pairs of points is O(D.N? anyway,
% so the improvement of selection is larger the smaller D is.
%
% In:
%   X: NxD matrix of N row D-dimensional vectors.
%   k: number of nearest neighbors.
%   method: find nearest neighbors with sorting ('sort') or selection ('mink').
%      Default: 'sort'.
% Out:
%   D2: N x k matrix of sorted square distances to the k nearest neighbors.
%   nn: N x k matrix of indices of the corresponding nearest neighbors.
%
% Any non-mandatory argument can be given the value [] to force it to take
% its default value.

% Copyright (c) 2013 by Max Vladymyrov and Miguel A. Carreira-Perpinan

function [D2,nn] = nnsqdist(X,k,method)

% ---------- Argument defaults ----------
if ~exist('method','var') || isempty(method) method='sort'; end
% ---------- End of "argument defaults" ----------

N = size(X,1); k = min(N-1,k);
% Block size (see below), select as large as your memory allows
mem = 1; B = floor((mem*1024^3)/(4*N*8)/2);	% This will fit in mem GB RAM

% Process by blocks to save memory
i1 = 1; i2 = min(N,B);
Xt = X'; X2 = X.^2; x2 = sum(X2,2)'; D2 = zeros(N,k); nn = D2;

while i1 <= N
  % This computes squared distances in a fast, vectorized way but can have
  % cancellation error for points closer than sqrt(eps).
  sd1 = max(bsxfun(@plus,sum(X2(i1:i2,:),2),...
                   bsxfun(@minus,x2,2*X(i1:i2,:)*Xt)),0);
  if method == 'sort'
    [sd1,ind] = sort(sd1,2);
  else
    [sd1,ind] = mink(sd1,k+1,2);
  end
  D2(i1:i2,:) = sd1(:,2:(k+1)); nn(i1:i2,:) = ind(:,2:(k+1));
  i1 = i1 + B; i2 = min(N,i1+B-1);
end
