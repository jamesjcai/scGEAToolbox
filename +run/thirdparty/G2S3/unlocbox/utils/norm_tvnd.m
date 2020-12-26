function y = norm_tvnd(u,type,weights)
%NORM_TVND N Dimentional TV norm
%   Usage:  norm_tvnd(x,weights)
%
%   Input parameters:
%         x     : Input data (N dimentional matrix)
%         type  : type ('isotropic' or 'anisotropic') (default 'isotropic')
%         weights: Weights
%   Output parameters:
%         sol   : Norm
%
%   Compute the N-dimentional TV norm of x
%
%   See also: norm_tv norm_tv3d
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/norm_tvnd.php

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


    if nargin < 2
        type = 'isotropic';
    end

    sz = size(u);
    dim = length(sz);
    
    if nargin<3
        weights = ones(dim,1);
    end


    temp = zeros(sz);

    if strcmp(type,'anisotropic')
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad)+ temp;
        end
        y = sum(temp(:));
    elseif strcmp(type,'isotropic')
        for d = 1:dim
            sz_temp = sz;
            sz_temp(d) = 1;
            tv(d).grad = weights(d)*cat(d,diff(u,1,d),zeros(sz_temp));
            temp = abs(tv(d).grad.^2)+ temp;
        end
        temp = sqrt(temp);
        y = sum(temp(:));
    else
        error('NORMTV_ND: unknown type.')
    end
end

