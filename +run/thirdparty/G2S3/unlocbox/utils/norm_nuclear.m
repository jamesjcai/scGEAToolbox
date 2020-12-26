function n = norm_nuclear(x)
%NORM_NUCLEAR - Nuclear norm of x
%   Usage: norm_nuclear(x) 
%  
%   Input parameters
%       x       : a matrix
%   Output parameters
%       n       : nuclear norm of x
%
%   See also: prox_nuclearnorm
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/norm_nuclear.php

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

% Author:  Vassilis Kalofolias, Nathanael Perraudin
% Date: February 2014
%

if issparse(x)
    n = sum(svds(x, min(size(x))));
else
    n = sum(svd(x, 'econ'));
end

end

