function dcor = distcorr(x, y)

% This function calculates the distance correlation between x and y.
% Reference: http://en.wikipedia.org/wiki/Distance_correlation
% Date: 18 Jan, 2013
% Author: Shen Liu (shen.liu@hotmail.com.au)

% Check if the sizes of the inputs match
if size(x, 1) ~= size(y, 1)
    error('Inputs must have the same number of rows')
end

% Delete rows containing unobserved values
% N = any([isnan(x) isnan(y)],2);
% x(N,:) = [];
% y(N,:) = [];

% Calculate doubly centered distance matrices for x and y


a = xpdist2(x);
A = a - mean(a) - mean(a, 2) + mean(a(:));
% mcol = mean(a);
% mrow = mean(a,2);
% ajbar = ones(size(mrow))*mcol;
% akbar = mrow*ones(size(mcol));
% abar = mean(mcol)*ones(size(a));
% A = a - ajbar - akbar + abar;

b = xpdist2(y);
B = b - mean(b) - mean(b, 2) + mean(b(:));

%mrow = mean(b,2);

% mcol = mean(b);
% mrow = mean(b,2);
% bjbar = ones(size(mrow))*mcol;
% bkbar = mrow*ones(size(mcol));
% bbar = mean(mcol)*ones(size(b));
% B = b - bjbar - bkbar + bbar;

% Calculate squared sample distance covariance and variances
dcov = sum(sum(A.*B)) / (size(x, 1)^2);

dvarx = sum(sum(A.*A)) / (size(x, 1)^2);
dvary = sum(sum(B.*B)) / (size(x, 1)^2);

% Calculate the distance correlation
dcor = sqrt(dcov/sqrt(dvarx*dvary));

end

function [D] = xpdist2(data)
d = dot(data, data, 2);
D = sqrt(d+d'-2*(data * data'));
end
