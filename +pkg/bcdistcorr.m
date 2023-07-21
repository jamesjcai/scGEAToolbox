function [bcR, p, T, df] = bcdistcorr(x, y)
% BCDISTCORR computes the bias corrected distance correlation
%
%   [BCR,P,T,DF] = BCDISTCORR(X,Y), where X (size n-by-p) and Y (size
%   n-by-q) are n random samples of observation.  The function returns the
%   bias corrected distance correlation BCR and the corresponding p-value
%   P, as well as the student t statistics T and its degree of freedom DF.
%   Note that the The t-test of independence is unbiased for every n ≥ 4
%   and any significance level.
%
%   This implementation is based on Székely, G. J., & Rizzo, M. L. (2013).
%   The distance correlation t-test of independence in high dimension.
%   Journal of Multivariate Analysis, 117, 193-213.
%
%   Date: 7.30.2016
%   Author: Po-He Tseng (pohetsn@gmail.com)

assert(rows(x)==rows(y));

n = rows(x);
X = Astar(x);    
Y = Astar(y);

XY = modified_distance_covariance(X, Y);
XX = modified_distance_covariance(X, X);
clear X;
YY = modified_distance_covariance(Y, Y);
clear Y;

bcR = XY/sqrt(XX*YY);
M = n*(n-3)/2;
T = sqrt(M-1) * bcR / sqrt(1-bcR^2);
df = M-1;
p = 1 - tcdf(T, df);

fprintf('bias-corrected R = %.3f, p-value=%.3f, T(%d)=%.4f\n',...
    bcR, p, df, T);

end

function XY = modified_distance_covariance(X, Y)
n = rows(X);
XY = sum(sum(bsxfun(@times, X, Y)))...
    - (n/(n-2))*diag(X)'*diag(Y);
end

function A = Astar(x)
d = pdist2(x,x);

n = rows(x);
m = mean(d);
M = mean(d(:));

% A = d - m'*ones(1,n);
A = bsxfun(@minus, d, bsxfun(@mtimes, m', ones(1,n)));
% A = A - ones(n,1)*m;
A = bsxfun(@minus, A, bsxfun(@mtimes, ones(n,1), m));
% A = A + M;
A = bsxfun(@plus, A, M); 

A = A - d/n;
A(1:n+1:end) = m - M;
A = (n/(n-1))*A;
end

function r = rows(x)
    r = size(x,1);
end