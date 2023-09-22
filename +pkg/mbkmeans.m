function [c, counts, idx] = mbkmeans(x, k, c, counts)

% [c1, counts1] = mbkmeans(subset1,numberofclusters);
% [c2, counts2] = mbkmeans(subset2,numberofclusters, c1, counts1); %start clustering using previously created clusters
% [c3, counts3] = mbkmeans(subset3,numberofclusters, c2, counts2);
% ...
% [cn, countsn, indices] = mbkmeans(subsetn,numberofclusters, c(n-1), counts(n-1));

% @MISC {164725,
%     TITLE = {How to apply distance-based clustering or dimensionality reduction for too many samples},
%     AUTHOR = {rcpinto (https://stats.stackexchange.com/users/49196/rcpinto)},
%     HOWPUBLISHED = {Cross Validated},
%     NOTE = {URL:https://stats.stackexchange.com/q/164725 (version: 2017-05-23)},
%     EPRINT = {https://stats.stackexchange.com/q/164725},
%     URL = {https://stats.stackexchange.com/q/164725}
% }

[N, D] = size(x);
if ~exist('c', 'var') || isempty(c)
    c = x(1:min([k, N]), :) + bsxfun(@times, randn(min([k, N]), D)*0.001, std(x));
    if N < k
        c(N+1:k, :) = bsxfun(@plus, mean(x), bsxfun(@times, randn(k-N, D), std(x)));
    end
end
if ~exist('counts', 'var') || isempty(counts)
    counts = zeros(k, 1);
end
idx = knnsearch(c, x, 'k', 1);
add = full(sparse(idx, 1, 1));
counts(idx) = counts(idx) + add(idx);
lr = 1 ./ counts(idx);
for i = 1:N
    c(idx(i), :) = (1 - lr(i)) * c(idx(i), :) + lr(i) * x(i, :);
end
end