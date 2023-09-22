function obj = filtadjc(obj, q)
if nargin < 2, q = 0.95; end
a = mean(maxk(abs(obj.A(:)), 10));
if a == 0, a = 1; end
obj.A = obj.A ./ a;
obj.A = obj.A .* (abs(obj.A) > quantile(abs(obj.A(:)), q));
if ~issparse(obj.A)
    obj.A = sparse(obj.A);
end
end