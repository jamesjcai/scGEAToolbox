function n0 = norm0(x,tol)
%NORM0 Compute the 0 norm of a vector
%   Usage: n0 = norm0(x);
%          n0 = norm0(x,tol);
%
%   Input parameters
%       x   : Vector or matrix
%       tol : Tolerance (default 1e-10)
%
%   Ouput parameters:
%       n0  : Zero norm
%
%   This function compute the zero norm of a vector or a matrix. More
%   explicitely, it count the number or non zero entries.
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/norm0.php

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
% Date  : 1 July 2014


if nargin<2
    tol = 1e-10;
end

n0 = sum(abs(x(:))>tol);

end
