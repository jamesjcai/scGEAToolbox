function [S, St] = sum_squareform(n, mask)
%SUM_SQUAREFORM sparse matrix that sums the squareform of a vector
%   Usage:  [S, St] = sum_squareform(n)
%           [S, St] = sum_squareform(n, mask)
%
%   Input parameters:
%         n:    size of matrix W
%         mask: if given, S only contain the columns indicated by the mask
%
%   Output parameters:
%         S:    matrix so that S*w = sum(W) for vector w = squareform(W)
%         St:   the adjoint of S
%
%   Creates sparse matrices S, St = S' so that
%       S*w = sum(W),       where w = squareform(W)
%
%   The mask is used for large scale computations where only a few
%   non-zeros in W are to be summed. It needs to be the same size as w,
%   n(n-1)/2 elements. See the example below for more details of usage.
%
%   Properties of S:
%    size(S) = [n, (n(n-1)/2)]     % if no mask is given.
%    size(S, 2) = nnz(w)           % if mask is given
%    norm(S)^2 = 2(n-1)
%    sum(S) = 2*ones(1, n*(n-1)/2)
%    sum(St) = sum(squareform(mask))   -- for full mask = (n-1)*ones(n,1)
%
%   Example:
%           % if mask is given, the resulting S are the ones we would get with the
%           % following operations (but memory efficiently):
%           [S, St] = sum_squareform(n);
%           [ind_i, ~, w] = find(mask(:));
%           % get rid of the columns of S corresponding to zeros in the mask
%           S = S(:, ind_i);
%           St = St(ind_i, :);
%
%   See also: squareform_sp
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/sum_squareform.php

% Copyright (C) 2012-2016 Nathanael Perraudin.
% This file is part of UNLOCBOX version 1.7.4
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%
% code author: Vassilis Kalofolias
% date: June 2015
% Testing: test_squareform




if nargin < 2
    mask_given = false;
else
    mask_given = true;
end


if mask_given
    %% more efficient than the following:
    % M = squareform(mask);
    % [I, J] = find(triu(M));
    
    if not(length(mask) == n*(n-1)/2)
        error('mask size has to be n(n-1)/2');
    end
    
    % TODO: make more efficient! (Maybe impossible)
    if iscolumn(mask)
        [ind_vec, ~] = find(mask);
    else
        [~, ind_vec] = find(mask);
    end
    
    % final number of columns is the nnz(mask)
    ncols = length(ind_vec);
    
    % indices inside the matrix
    I = zeros(ncols, 1);
    J = zeros(ncols, 1);
    
    curr_row = 1;
    offset = 0;
    % length of current row of matrix, counting from after the diagonal
    len = n - 1;
    for ii = 1 : ncols
        ind_vec_i = ind_vec(ii);
        % if we change row, the vector index is bigger by at least the
        % length of the line + the offset.
        while(ind_vec_i > (len + offset))
            offset = offset + len;
            len = len - 1;
            curr_row = curr_row + 1;
        end
        I(ii) = curr_row;
        J(ii) = ind_vec_i - offset + (n-len);
    end
    
else
    %% more efficient than the following:
    % W = ones(n) - eye(n);
    % [I, J] = find(tril(W));
    
    % number of columns is the length of w given size of W
    ncols = (n-1)*(n)/2;
    
    I = zeros(ncols, 1);
    J = zeros(ncols, 1);
    
    % offset
    k = 1;
    for i = 2 : n
        I(k: k + (n-i)) = i : n;
        k = k + (n-i+1);
    end
    k = 1;
    for i = 2 : n
        J(k: k + (n-i)) = i-1;
        k = k + (n-i+1);
    end
end

St = sparse([1:ncols, 1:ncols], [I, J], 1, ncols, n);
S = St';

end
