function [ sol]  = prox_add_2norm( x,gamma,param )
%PROX_ADD_2NORM Proximal operator with an additional quadratic term
%   Usage:   sol = prox_add_2norm(x, gamma, param);
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         f     : Function
%   Output parameters:
%         sol   : Solution.
%         infos : Structure summarizing informations at convergence
%
%   PROX_ADD_2NORM( x,gamma,param ) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + 0.5*||y - z||_2^2 + gamma * f(z)
%
%   This problem can be solved because we have the nice relationship
%
%      0.5*||x - z||_2^2 + 0.5*||y - z||_2^2 = || (x+y)/2 - z||_2^2 
%                                              + 0.25 *||y - x||_2^2
%
%   This can be used to reduce the number of functionals and the solution is 
%
%      sol = prox_{gamma/2 * f} ((x+y)/2)
%
%   param is a Matlab structure containing the following fields:
%   
%    param.y : a vector of the same size as x
%
%    param.f : a structure containing the function f
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/prox_add_2norm.php

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
% Date: 1 February 2014
%



% Optional input arguments
if nargin<3,     
    error('Two few input arguments!');
end

if ~isfield(param, 'y')
   error('Please specify param.y')
end
if ~isfield(param, 'f')
   error('Please specify param.f')
end

sol = param.f.prox((x+param.y)/2,0.5*gamma);

end


