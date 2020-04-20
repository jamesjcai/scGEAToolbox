function W = zero_diag(W)
%ZERO_DIAG sets the diagonal of a matrix to 0
%   Usage: B = zero_diag(A);
%
%   Input parameters:
%       A   : input matrix
%   Output parameters:
%       B   : output with zero diagonal
%
%   Works also for non-square matrices
%
%   See also: squareform_sp
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/zero_diag.php

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

% code author: Vassilis Kalofolias
% date: March 2016


[m, n] = size(W);

% if not(m==n)
%     warning('non square matrix given!')
% end

n_zeros = min(m, n);

W(1: (m+1): (m+1) * n_zeros) = 0;



