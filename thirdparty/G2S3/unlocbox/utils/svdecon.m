function [U,S,V] = svdecon(X)
%SVDECON Fast svds when n<<m
%   Usage: [U,S,V] = svdecon(X);
%
%   Input parameters:
%         X     : Input data (n x m)
%   Output parameters:
%         U     : Left singular vectors
%         S     : Singular values
%         U     : Right signular vectors
%
%   This function is an acceleration of svd. It is particularly efficient
%   when n<<m
%
%   See also: svdsecon
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/svdecon.php

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

    [m,n] = size(X);

    if  m <= n
        X2 = X*X';
        [U,E] = eig(X2);

        [e,ix] = sort(abs(diag(E)),'descend');
        U = U(:,ix);    

        V = X'*U;
        s = sqrt(e);
        
        V = bsxfun(@times, V, 1./s');
        S = diag(s);
    else
        X2 = X'*X; 
        [V,E] = eig(X2);

        [e,ix] = sort(abs(diag(E)),'descend');
        V = V(:,ix);    
        %% TODO: compute U = X*(V/S) and NOT U = (X*V)/S
        U = X*V; 
        s = sqrt(e);
        U = bsxfun(@times, U, 1./s');
        S = diag(s);
    end
end
