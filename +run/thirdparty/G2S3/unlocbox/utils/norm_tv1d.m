function y = norm_tv1d(I,w)
%NORM_TV1D 1 Dimentional TV norm
%   Usage:  y = norm_tv1d(x)
%           y = norm_tv1d(x,w)
%
%   Input parameters:
%         I     : Input data 
%         w    : Weights
%   Output parameters:
%         y     : Norm
%
%   Compute the 1-dimentional TV norm of I. If the input I is a matrix.
%   This function will compute the norm of all line and return a vector of
%   norms.
%
%   See also: norm_tv norm_tv3d
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/norm_tv1d.php

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

if nargin>1
    dx = gradient_op1d(I,w);
else
    dx = gradient_op(I);
end    

%y = sum(temp(:));
y = sum(abs(dx),1)';

end

