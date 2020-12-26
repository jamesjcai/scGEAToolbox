function [sz] = soft_threshold(z,T)
%SOFT_THRESHOLD soft thresholding
%   Usage:  sz = soft_threshold(z,T);
%
%   Input parameters:
%         z     : Input signal
%         T     : Threshold
%                 if T is a vector, then thresholding is applied component-wise
%
%   Output parameters:
%         sz    : Soft thresholded signal
%   
%   This function soft thresholds z by T. It can handle complex input z.
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/soft_threshold.php

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

% Nathanael Perraudin
% Date: 14 May 2013


size_z = size(z);

if all(T==0)
    sz = z;  %identity
elseif any(T<0)
    error('Threshold value(s) cannot be negative!')
elseif isscalar(T)  %for scalar threshold it is faster to compute it like this
    % handle the size
    z=z(:);
    
    % This soft thresholding function only supports real signal
    % sz = sign(z).*max(abs(z)-T, 0);
    
    % This soft thresholding function supports complex numbers
    sz = max(abs(z)-T,0)./(max(abs(z)-T,0)+T).*z;
    
else  %for vector threshold(s) it is faster to compute it like this
    % handle the size
    z=z(:);
    T=T(:);
    
    % This soft thresholding function supports complex numbers
    % sz = max(abs(z)-T,0)./(max(abs(z)-T,0)+T).*z;
    aux = max(abs(z)-T,0);
    sz = aux./(aux+T).*z;

end

% Handle the size
sz = reshape(sz,size_z);

