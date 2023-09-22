function y = sparse_nanmean(X, dim)
if nargin < 2, dim = 1; end
if dim ~= 1, X = X'; end
Z = isnan(X) | (X == 0);
z = sum(~Z);
y = sum(X, 'omitnan');
y = y ./ z;
if dim ~= 1, y = y'; end
end
