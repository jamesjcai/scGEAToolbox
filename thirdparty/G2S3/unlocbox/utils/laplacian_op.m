function [ L ] = laplacian_op( I )
%LAPLACIAN_OP 2 dimentional Laplacian
%   Usage:  [I] = laplacian_op( I );
%
%   Input parameters:
%         I     : Input image 
%
%   Output parameters:
%         I     : Laplacian
%
%   Compute the sum of the laplacian along x and y. This operator is
%   self-adjoint.
%
%           L = I_xx + I_yy
%
%   See also: laplacianx_op laplaciany_op div_op gradient_op
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/laplacian_op.php

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

% Author: Nathnaael Perraudin
% Date  : 13 September 2013


L=laplacianx_op(I)+laplaciany_op(I);

end


