%NANGEOMEAN Geometric mean, ignoring NaNs.
%   M = NANGEOMEAN(X) returns the geometric mean of X, treating NaNs as
%   missing values. If X is a vector, M is the n-th root of the product of
%   the non-NaN elements. If X is a matrix, M is a vector of the geometric
%   mean of each column. For N-D arrays, NANGEOMEAN operates along the
%   first non-singleton dimension.
%   
%   M = NANGEOMEAN(X, DIM) takes the geometric mean along dimension DIM.
%   
%   See also NANMEAN, GEOMEAN.

function m = nangeomean(x, dim)

% Default inputs
if nargin<2
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Input sanitisation
if any(x(:) < 0)
    error(message('nangeomean:BadData'))
end

% Find NaNs and set them to one
nans = isnan(x);
x(nans) = 1;

% Count up non-NaNs.
n = sum(~nans,dim);
n(n==0) = NaN; % prevent divideByZero warnings

% Geometric mean is equivalent to the exponent of the mean in logspace
m = exp(sum(log(x),dim)./n);
% Note:
%   Since log(1)==0, the NaNs do not contribute to the sum.
%   Since n is the count of only non-NaN elements, the denominator is
%   correct.

end