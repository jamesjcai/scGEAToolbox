function I = div_op1d(dx, wx)
%DIV_OP1D Divergence operator in 1 dimention
%   Usage:  I = div_op1d(dx)
%           I = div_op1d(dx, wx)
%
%   Input parameters:
%         dx    : Gradient along x
%         wx    : Weights along x
%
%   Output parameters:
%         I     : Output divergence vector 
%
%   Compute the 1-dimentional divergence of a vector. If a matrix is given,
%   it will compute the divergence of all vectors in the matrix.
%
%   Warning this function compute the divergence operator defined as minus
%   the adjoint of the gradient
%
%           div  = - grad'
%
%   See also: gradient_op div_op3d div_op1d laplacian_op
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/div_op1d.php

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
end

I = [dx(1, :) ; dx(2:end-1, :)-dx(1:end-2, :) ; -dx(end-1, :)];

end

