function result_data = general_sset_intersection(left, right, result_row, result_col, result_data, mix_weight)
%GENERAL_SSET_INTERSECTION Calculate the values in a sparse matrix
% resulting from a simplicial set intersection between "left" and "right".
%
% result_data = GENERAL_SSET_INTERSECTION(left, right, result_row, result_col, result_data, mix_weight)
% 
% Parameters
% ----------
% left: sparse matrix
%     The first input fuzzy simplicial set.
% 
% right: sparse matrix
%     The second input fuzzy simplicial set.
%
% result_row: array of size (n_entries, 1)
%     The rows of non-zero entries in the result.
%
% result_col: array of size (n_entries, 1)
%     The columns of non-zero entries in the result.
%
% result_data: array of size (n_entries, 1)
%     The values resulting from the sum of "left" and "right".
% 
% mix_weight: double (optional, default 0.5)
%     The weight assigned to "left" when performing the
%     intersection.
% 
% Returns
% -------
% result_data: sparse matrix
%     The non-zero entries of the intersected fuzzy simplicial set.

%   AUTHORSHIP
%   Math Lead & Primary Developer:  Connor Meehan <connor.gw.meehan@gmail.com>
%   Secondary Developer: Stephen Meehan <swmeehan@stanford.edu>
%   Bioinformatics Lead:  Wayne Moore <wmoore@stanford.edu>
%   Copyright (c) 2022 The Board of Trustees of the Leland Stanford Junior University; Herzenberg Lab
%   License: BSD 3 clause

    if nargin < 6
        mix_weight = 0.5;
    end
    
    left_data = nonzeros(left);
    right_data = nonzeros(right);

    if isempty(left_data)
        left_min = 1e-8;
    else
        left_min = max(min(left_data)/2, 1e-8);
    end
    if isempty(right_data)
        right_min = 1e-8;
    else
        right_min = max(min(right_data)/2, 1e-8);
    end

    left_vals = full(left(sub2ind(size(left), result_row, result_col)));
    left_vals = left_vals + left_min*(left_vals == 0);
    right_vals = full(right(sub2ind(size(right), result_row, result_col)));
    right_vals = right_vals + right_min*(right_vals == 0);
    
    change_val = left_vals > left_min | right_vals > right_min;
    if mix_weight < 0.5
        result_data(change_val) = left_vals.*right_vals.^(mix_weight/(1 - mix_weight));
    else
        result_data(change_val) = right_vals.*left_vals.^(mix_weight/(1 - mix_weight));
    end
end
    