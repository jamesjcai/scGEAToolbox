function n21 = norm_l21(x, g_d,g_t, w2,w1)
%NORM_L21 L21 mixed norm
%   Usage:  n21 = norm_l21(x);
%           n21 = norm_l21(x, g_d,g_t);
%           n21 = norm_l21(x, g_d,g_t, w2,w1);
%
%   Input parameters:
%         x     : Input data 
%         g_d   : group vector 1
%         g_t   : group vector 2
%         w2    : weights for the two norm (default 1)
%         w1    : weights for the one norm (default 1)
%   Output parameters:
%         y     : Norm
%
%   NORM_L21(x, g_d,g_t, w2,w1) returns the norm L21 of x. If x is a
%   matrix the 2 norm will be computed as follow:
%
%       n21 = || x ||_21 = sum_j ( sum_i |x(i,j)|^2 )^(1/2) 
%
%   In this case, all other argument are not necessary.
%
%   'norm_l21(x)' with x a row vector is equivalent to norm(x,1) and
%   'norm_l21(x)' with x a line vector is equivalent to norm(x)
%
%   For fancy group, please provide the groups vectors.
%
%   g_d, g_t are the group vectors. g_d contain the indices of the
%   element to be group and g_t the size of different groups.
%       
%   Example: 
%                x=[x1 x2 x3 x4 x5 x6] 
%                Group 1: [x1 x2 x4 x5] 
%                Group 2: [x3 x6]
%
%   Leads to 
%           
%               => g_d=[1 2 4 5 3 6] and g_t=[4 2]
%               Or this is also possible
%               => g_d=[4 5 3 6 1 2] and g_t=[2 4]   
%
%   This function works also for overlapping groups.
%
%   See also: norm_linf1 norm_tv
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/norm_l21.php

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
% Date: October 2011
% Testing: test_mixed_sparsity


% Optional input arguments
if nargin<2, g_d = 1:numel(x); end
if nargin<3, 
    if numel(x) == size(x,1)*size(x,2); % matrix case
        g_t = size(x,2)*ones(1,size(x,1)); 
    else
        g_t = ones(1,numel(x)); 
    end
end
if nargin<5, w1=ones(size(g_t,2),size(g_t,1)); end
if nargin<4, w2=ones(numel(x),size(g_t,1)); end


% overlapping groups
if size(g_d,1)>1
    n21=0;
    for ii=1:size(g_d,1);
        n21 = n21 + norm_l21(x,g_d(ii,:),g_t(ii,:), ...
            w2(:,ii),w1(:,ii));
    end
else % non overlapping groups


    l=length(g_t);

    % Compute the norm
    if max(g_t)==min(g_t) % group of the same size
        X = transpose(x);
        X = X(g_d);
        X = transpose(reshape(X,numel(x)/l,l));

        W2 = transpose(reshape(w2(g_d),numel(x)/l,l));
        normX2 = sqrt(sum((abs(W2.*X)).^2,2));
        n21 = sum(w1.*normX2);
    else % group of different size

        n21 = 0;
        indice = 0;
        X = x(:);
        X = X(g_d);
        W2 = w2;
        W2 = W2(g_d);
        for i=1:l

            n21 = n21 + w1(i) * norm(W2(indice+1:indice+g_t(i)) ...
                         .*X(indice+1:indice+g_t(i)));
            indice = indice+g_t(i);
        end

    end
    
end

end

