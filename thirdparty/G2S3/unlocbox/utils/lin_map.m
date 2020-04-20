% LIN_MAP   Map linearly from a given range to another
%
% USAGE:
% Y = lin_map(X, lims_out, lims_in)
% 
% X can be of any size. Elements of lims_in or lims_out don't have to be in
% an ascending order.
%
% if lims_in is not specified, the minimum and maximum values are used. 
%
% Example: 
%     
% x = cos((1:50)/3) + .05*randn(1, 50);
% y = lin_map(x, [4, 2]);
% figure; plot(x, 'o'); hold all; plot(y, '*');
%
%
%
%
%code author: Vassilis Kalofolias
%date: Feb 2014
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/lin_map.php

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





function [Y] = lin_map(X, lims_out, lims_in)

if nargin < 3
    lims_in = [min(X(:)), max(X(:))];
end
    

a = lims_in(1);
b = lims_in(2);
c = lims_out(1);
d = lims_out(2);


Y = zeros(size(X), class(X));


Y(:) = (X(:)-a) * ((d-c)/(b-a)) + c;




    

    
    
    
