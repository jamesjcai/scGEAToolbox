function [ Ly ] = laplaciany_op( I )
%LAPLACIANY_OP dimentional Laplacian
%   Usage:  [Ly] = laplaciany_op( I );
%
%   Input parameters:
%         I     : Input image 
%
%   Output parameters:
%         Ly    : Laplacian along y
%
%   Compute the sum of the laplacian along y. This operator is
%   self-adjoint.
%
%           Ly = I_yy
%
%   See also: laplacian_op laplacianx_op div_op gradient_op
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/laplaciany_op.php

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


dy = [I(:, 2:end,:)-I(:, 1:end-1,:) , zeros(size(I, 1), 1,size(I, 3))];

Ly =  [dy(:, 1,:) , dy(:, 2:end-1,:)-dy(:, 1:end-2,:) , -dy(:, end-1,:)];

end

