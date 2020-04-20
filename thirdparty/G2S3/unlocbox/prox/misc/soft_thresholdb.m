function [sz] = soft_thresholdb(z,T)
%SOFT_THRESHOLDB soft thresholding for mixed sparsity
%   Usage:  soft_threshold(z,T)
%
%   Input parameters:
%         z     : Input signal.
%         T     : Threshold.
%   Output parameters:
%         sz    : Soft thresholded signal.
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/misc/soft_thresholdb.php

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

    % handle the size
    size_z=size(z);
    z=z(:);
    T=T(:);
    
    % Precaution on T
    T(T==inf)=realmax;
    

    % This soft thresholding function support complex number
    sz = max(abs(z)-T.*abs(z),0)./(max(abs(z)-T.*abs(z),0)+T.*abs(z)+double(abs(z)<eps)).*z;
    

    % Handle the size
    sz=reshape(sz,size_z);
end

