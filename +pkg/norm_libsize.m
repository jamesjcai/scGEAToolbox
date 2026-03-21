function X = norm_libsize(X, factor)
% norm_libsize  Normalize count matrix by library size and scale, returning normalized data.
%   Xnorm = norm_libsize(X) normalizes each column to counts per 10,000,
%   preserving double type. Optionally log1p(Xnorm) can be applied externally.
%
%   Xnorm = norm_libsize(X, factor) uses custom scaling factor; default = 1e4.

% ensure double
X = double(X);
libSize = sum(X, 1, 'omitnan');  % 1×n vector

if nargin < 2 || isempty(factor)
    factor = 1e4;
end
libSize(libSize == 0) = NaN;
X = (X ./ libSize) * factor;
end
