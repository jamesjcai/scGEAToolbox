function set_seed(my_seed)
%SET_SEED sets the seed of the default random random generator
%   Usage:  set_seed(my_seed)
%           set_seed()
%
%   Input parameters:
%         my_seed  :  new_seed
%   
%   Set the seed of the default random random generator
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/set_seed.php

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

% code author: Vassilis Kalofolias
% date: August 2013


%% Fix the random stream for debugging reasons
if nargin < 1
    my_seed = 0;
end

if verLessThan('matlab', '7.12.0')  % release 2011a has "rng"
    rand('twister', my_seed);
else
    rng(my_seed, 'twister');
    %rng('default');
end

end




