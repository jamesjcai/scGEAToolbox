function [row,col,val]=i_maxij(matrix,k)
    % Flatten the matrix and find the top 10 maximum values
    [sorted_values, sorted_indices] = sort(matrix(:), 'descend');
    top_10_indices = sorted_indices(1:k);
    % Convert linear indices to row and column subscripts
    [row, col] = ind2sub(size(matrix), top_10_indices);
    if nargout > 2
        val = sorted_values(1:k);
    end
end
