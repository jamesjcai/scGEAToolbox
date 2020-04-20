function I = div_op(dx, dy, wx, wy)
%DIV_OP Divergence operator in 2 dimensions
%   Usage:  I = div_op(dx, dy)
%           I = div_op(dx, dy, wx, wy)
%
%   Input parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%         wx    : Weights along x
%         wy    : Weights along y
%
%   Output parameters:
%         I     : Output divergence image 
%
%   Compute the 2-dimensional divergence of an image. If a cube is given,
%   it will compute the divergence of all images in the cube.
%
%   Warning: computes the divergence operator defined as minus the adjoint
%   of the gradient 
%
%           div  = - grad'
%
%   See also: gradient_op div_op3d div_op1d laplacian_op prox_tv
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/div_op.php

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

if nargin > 2
    dx = dx .* conj(wx);
    dy = dy .* conj(wy);
end

I = [dx(1, :,:) ; ...
    dx(2:end-1, :,:)-dx(1:end-2, :,:) ;...
    -dx(end-1, :,:)];
I = I + [dy(:, 1,:) ,...
    dy(:, 2:end-1,:)-dy(:, 1:end-2,:) ,...
    -dy(:, end-1,:)];

end

