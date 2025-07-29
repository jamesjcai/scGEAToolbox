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
% Xnorm = bsxfun(@rdivide, X, libSize) * factor;
X = (X ./ libSize) * factor;
end


%{
function [X] = norm_libsize(X, factorn)
%Library size normalization

%aka: namely cell depth normalization to the mean cell depth, followed by log1p (log1pPF)
% https://www.biorxiv.org/content/10.1101/2022.05.06.490859v1.full
if ~isa(X, 'double')
    X = double(X);
end
lbsz = sum(X, 'omitnan');
if nargin < 2 || isempty(factorn)
    % factorn = mean(lbsz, 'omitnan');
    factorn = 1e4;
end
X = (X ./ lbsz) * factorn;
end

%}