function [A] = e_filtadjc(A, q, makesparse)
% Filter adjacnecy matrix A with cutoff Q
%
% A - adjacency matrix

if nargin < 3, makesparse = true; end
if nargin < 2, q = 0.95; end
dim = size(A);
if numel(dim) == 2
    a = mean(maxk(abs(A(:)), 10)); % top 10 average
    if a > 0
        A = A ./ a;
        A = A .* (abs(A) > i_quantile(abs(nonzeros(A)), q));
        if ~issparse(A) && makesparse
            A = sparse(A);
        end
    end
elseif numel(dim) == 3
    for k = 1:dim(3)
        A(:, :, k) = ten.e_filtadjc(A(:, :, k), q);
    end
end

end


function t = i_quantile(v, q)
% quantile(v, q) for a column v and scalar q, without sorting all of v.
%
% quantile sorts every value, 0.15 s for a 3000-gene network and so a
% sixth of the time scTenifoldNet spends on each one. Only two order
% statistics matter: with n values, quantile interpolates linearly between
% the k-th and (k+1)-th smallest, k = floor(n*q + 0.5). So bracket them
% from a fixed-stride sample (no random draws), keep the values inside the
% bracket, and sort only those. If the bracket misses either order
% statistic, which a sample can, fall back to quantile itself; the result
% is quantile's either way, bit for bit.
n = numel(v);
r = n*q + 0.5;
k = floor(r);
minSample = 1e5;
if n < 10*minSample || k < 1 || k >= n || any(isnan(v))
    t = quantile(v, q);
    return
end
stride = floor(n/minSample);
s = sort(v(1:stride:end));
ns = numel(s);
margin = 0.005;   % about 7 standard errors of a quantile from 1e5 values
lo = s(max(1, floor(ns*(q - margin))));
hi = s(min(ns, ceil(ns*(q + margin))));
below = nnz(v < lo);
window = sort(v(v >= lo & v <= hi));
kw = k - below;
if kw < 1 || kw + 1 > numel(window)
    t = quantile(v, q);
    return
end
t = window(kw) + (r - k)*(window(kw + 1) - window(kw));
end
