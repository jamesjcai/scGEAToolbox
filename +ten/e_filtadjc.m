function [A] = e_filtadjc(A, q, makesparse)
%Filter adjacnecy matrix A with cutoff Q
%
% A - adjacency matrix

if nargin < 3, makesparse = true; end
if nargin < 2, q = 0.95; end
dim = size(A);
if numel(dim) == 2
    a = mean(maxk(abs(A(:)), 10)); % top 10 average
    % a=max(abs(A(:)));
    if a > 0
        A = A ./ a;
        A = A .* (abs(A) > quantile(abs(nonzeros(A)), q));
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
