function I = div_op3d(dx, dy, dz, wx, wy, wz)
%DIV_OP3D Divergence operator in 3 dimentions
%   Usage:  I = div_op3d(dx, dy, dz)
%           I = div_op3d(dx, dy, dz, wx, wy, wz)
%
%   Input parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%         dz    : Gradient along z
%         wx    : Weights along x
%         wy    : Weights along y
%         wz    : Weights along z
%
%   Output parameters:
%         I     : Output image
%
%   Compute the 3-dimentional divergence of a 3D-image. If a 4 dimentional
%   signal is given, it will compute the divergence of all cubes in the
%   4 diementionals signal.  
%
%   Warning this function compute the divergence operator defined as minus
%   the adjoint of the gradient
%
%           div  = - grad'
%
%   See also: gradient_op div_op div_op1d laplacian_op
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/div_op3d.php

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

if nargin > 3
    dx = dx .* conj(wx);
    dy = dy .* conj(wy);
    dz = dz .* conj(wz);
end

I = [dx(1, :, :,:) ; dx(2:end-1, :, :,:) - ...
    dx(1:end-2, :, :,:) ; -dx(end-1, :, :,:)];
I = I + [dy(:, 1, :,:) , dy(:, 2:end-1, :,:) - ...
    dy(:, 1:end-2, :,:) , -dy(:, end-1, :,:)];
I = I + cat(3, dz(:, :, 1,:) , dz(:, :, 2:end-1,:) - ...
    dz(:, :, 1:end-2,:) , -dz(:, :, end-1,:));
end

