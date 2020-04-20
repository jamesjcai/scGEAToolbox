function [ stop ]=test_gamma(gamma)
%TEST_GAMMA test if gamma is correct
%   Usage:  stop = test_gamma(gamma)
%           test_gamma(gamma)
%
%   Input parameters:
%         gamma : number
%   Output parameters:
%         stop  : boolean
%
%   This function test is gamma is stricly positive
%   
%   If gamma is negativ, this function return an error. If gamma is zero
%   this function, set stop to 1.
%
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/misc/test_gamma.php

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

% Author:  Nathanael Perraudin
% Date: February 2012
%


if gamma<0
    error('gamma can not be negativ!');
% elseif (gamma==0) && warning
%     gamma=gamma+eps;
%     fprintf(' WARNING!!! gamma is 0. We add eps to gamma to keep going...\n');
% else
%   % gamma = gamma; 
end
    
if gamma==0
    stop = 1;
else
    stop = 0;
end

end

