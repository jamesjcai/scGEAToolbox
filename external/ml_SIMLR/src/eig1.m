function [eigvec, eigval, eigval_full] = eig1(A, c, isMax, isSym)

if nargin < 2
    c = size(A, 1);
    isMax = 1;
    isSym = 1;
elseif c > size(A, 1)
    c = size(A, 1);
end

if nargin < 3
    isMax = 1;
    isSym = 1;
end

if nargin < 4
    isSym = 1;
end

if isSym == 1
    A = max(A, A');
end
[v, d] = eig(A);
d = diag(d);
%d = real(d);
if isMax == 0
    [~, idx] = sort(d);
else
    [~, idx] = sort(d, 'descend');
end

idx1 = idx(1:c);
eigval = d(idx1);
%eigvec = real(v(:,idx1)*diag(sqrt(1-eigval)));
eigvec = real(v(:, idx1));
eigval_full = d(idx);