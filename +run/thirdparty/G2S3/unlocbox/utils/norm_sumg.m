function n = norm_sumg(x, G, w)
%NORM_SUMG 2 Dimentional TV norm
%   Usage:  y = norm_sumg(x, G);
%           y = norm_sumg(x, G, w);
%
%   Input parameters:
%         x     : Input data (vector)
%         G     : The structure array of norm operator: 
%         w     : Weights (default 1)
%   Output parameters:
%         n     : Norm
%
%   n = NORM_SUMG(x, G, w) returns the sum of the norm x given in
%   the structure array G. The norm can be weighted using the parameter
%   weights.
%
%   See also: prox_sumg
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/norm_sumg.php

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
% Date: October 2011
%


% Optional input arguments
if nargin<2, error('No input functions!'); end
if nargin<3, w=ones(length(G),1); end


% Compute the norm
n=0;

for ii=1:length(G)
    n=n+w(ii)*G{ii}.eval(x);
end

end

