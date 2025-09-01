%{
% https://github.com/weihuayi/numopde/blob/master/sparse.md
X = [0 0 0 10; 21 0 33 0; 0 0 3 0; 12 1 0 4];
X = sparse(X);

% pkg.e_writeh5(X,["aaa"],'ccc.h5');

[indptr, indices, data] = convert_sparse_to_indptr(X);
% indices = indices'-1    % J


%k = 1;
%data(indptr(1):indptr(k+1))
Y = zeros(size(X));

for k = 1:length(indptr)-1
    ix = indptr(k) + 1:indptr(k+1);
    ix = ix - 1;
    y = indices(ix);
    Y(y, k) = data(ix);
end

size(Y)
isequal(X,Y)




function [indptr, indices, data] = convert_sparse_to_indptr(X)

% Check if X is sparse
if ~issparse(X)
    error('Input matrix X must be a sparse matrix');
end

% Get matrix dimensions
[~, n] = size(X);

% Initialize indptr with 1 and n+1
indptr = [1, n+1];

% Find non-zero elements and their indices
[row, col] = find(X);

% Sort by columns for efficient construction
[~, sort_idx] = sort(col);
row = row(sort_idx);
col = col(sort_idx);

% Accumulate column counts for indptr
for i = 1:n
    indptr(i+1) = indptr(i) + sum(col == i);
end

% Assign indices and data
indices = row;
% data = full(X(row, col));  % Extract non-zero values
data = nonzeros(X);

end
%}