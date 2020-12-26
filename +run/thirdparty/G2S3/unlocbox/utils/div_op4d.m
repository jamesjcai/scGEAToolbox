function I = div_op4d(dx, dy, dz, dt, wx, wy, wz, wt)
%DIV_OP4D Divergence operator in 4 dimentions
%   Usage:  I = div_op4d(dx, dy, dz, dt)
%           I = div_op4d(dx, dy, dz, dt, wx, wy, wz, wt)
%
%   Input parameters:
%         dx    : Gradient along x
%         dy    : Gradient along y
%         dz    : Gradient along z
%         dt    : Gradient along t
%         wx    : Weights along x
%         wy    : Weights along y
%         wz    : Weights along z
%         wt    : Weights along t
%
%   Output parameters:
%         I     : Output image
%
%   Compute the 4-dimentional divergence of a 4D-image. If a 5 dimentional
%   signal is given, it will compute the divergence of all 4 dimentional
%   cubes in the 5 diementionals signal.  
%
%   Warning this function compute the divergence operator defined as minus
%   the adjoint of the gradient
%
%           div  = - grad'
%
%   See also: gradient_op4d div_op div_op1d div_op3d
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/div_op4d.php

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

if nargin > 4
    dx = dx .* conj(wx);
    dy = dy .* conj(wy);
    dz = dz .* conj(wz);
    dt = dt .* conj(wt);
end

I = [dx(1, :, :, :, :) ; dx(2:end-1, :, :, :, :) - ...
    dx(1:end-2, :, :, :, :) ; -dx(end-1, :, :, :, :)];
I = I + [dy(:, 1, :, :, :) , dy(:, 2:end-1, :, :, :) - ...
    dy(:, 1:end-2, :, :, :) , -dy(:, end-1, :, :, :)];
I = I + cat(3, dz(:, :, 1, :, :) , dz(:, :, 2:end-1, :, :) - ...
    dz(:, :, 1:end-2, :, :) , -dz(:, :, end-1,:, :));
I = I + cat(4, dt(:, :, :, 1, :) , dt(:, :, :, 2:end-1, :) - ...
    dt(:, :, :, 1:end-2, :) , -dt(:, :, :, end-1, :));
end

