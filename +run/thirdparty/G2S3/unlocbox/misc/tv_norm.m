function y = tv_norm(I,wx,wy)
%TV_NORM 2 Dimentional TV norm
%   Usage:  y = tv_norm(x)
%
%   Input parameters:
%         I     : Input data 
%         wx    : Weights along x
%         wy    : Weights along y
%   Output parameters:
%         y     : Norm
%
%   Compute the 2-dimentional TV norm of I. If the input I is a cube. This
%   function will compute the norm of all image and return a vector of
%   norms.
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/misc/tv_norm.php

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

warning('This function will be deleted in future release of matlab, use norm_tv instead')

if nargin>1
    [dx, dy] = gradient_op(I,wx, wy);
else
    [dx, dy] = gradient_op(I);
end    
temp = sqrt(abs(dx).^2 + abs(dy).^2);

%y = sum(temp(:));
y = reshape(sum(sum(temp,1),2),[],1);

end

