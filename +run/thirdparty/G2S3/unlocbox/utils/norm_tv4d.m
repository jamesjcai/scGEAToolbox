function y = norm_tv4d(u,wx, wy, wz, wt)
%NORM_TV4D 4 Dimentional TV norm
%   Usage:  y = norm_tv4d(x)
%           y = norm_tv4d(x, wx, wy, wz, wt )
%
%   Input parameters:
%         x     : Input data (3 dimentional matrix)
%         wx    : Weights along x
%         wy    : Weights along y
%         wz    : Weights along z
%         wt    : Weights along t
%
%   Output parameters:
%         y   : Norm
%
%   Compute the 4-dimentional TV norm of x. If the input I is a 5
%   dimentional signal. This function will compute the norm of all 4
%   dimentional cubes and return a vector of norms.
%
%   See also: norm_tv norm_tvnd norm_tv3d
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/norm_tv4d.php

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
% Date:   24 April 2014

if nargin>1
    [dx, dy, dz, dt] = gradient_op4d(u,wx, wy, wz, wt);
else
    [dx, dy, dz, dt] = gradient_op4d(u);
end
    
temp = sqrt(abs(dx).^2 + abs(dy).^2 + abs(dz).^2+ abs(dt).^2);
% y = sum(temp(:));

% This allows to return a vector of norms
y = reshape(sum(sum(sum(sum(temp,1),2),3),4),[],1);


end

