function [dx, dy] = gradient_op(I, wx, wy)
%GRADIENT_OP 2 Dimensional gradient operator
%   Usage:  [dx, dy] = gradient_op(I)
%           [dx, dy] = gradient_op(I, wx, wy)
%
%   Input parameters:
%         I     : Input data 
%         wx    : Weights along x
%         wy    : Weights along y
%
%   Output parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%
%   Compute the 2-dimensional gradient of I. If the input I is a cube. This
%   function will compute the gradient of all image and return two cubes.
%
%   See also: gradient_op3d gradient_op1d div_op laplacian_op
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/gradient_op.php

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

% Author: Nathanael Perraudin
% Date:   1 February 2014

dx = [I(2:end, :,:)-I(1:end-1, :,:) ; zeros(1, size(I, 2),size(I, 3))];
dy = [I(:, 2:end,:)-I(:, 1:end-1,:) , zeros(size(I, 1), 1,size(I, 3))];

if nargin>1
    dx = dx .* wx;
    dy = dy .* wy;
end

end

