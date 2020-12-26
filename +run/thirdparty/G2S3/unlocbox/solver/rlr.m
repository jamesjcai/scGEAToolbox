function [sol, info] = rlr(x_0,f,A,At, param)
%RLR Regularized Linear Regression 
%   Usage: sol = rlr(x_0,f,A,At, param)
%          sol = rlr(x_0,f,A,At)
%          [sol, info] = rlr(..,)
%
%   Input parameters:
%         x_0   : Starting point of the algorithm
%         f     : Function to minimize
%         A     : Operator
%         At    : Adjoint operator
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%
%   This function solve minimization problem using forward-backward splitting
%
%   sol = RLR(x_0,f,A,At, param) solves:
%
%      sol = argmin ||x_0-Ax||_2^2 + f2(x)      for x belong to R^N
%
%   where x is the variable. 
%
%    x_0 is the starting point.
%
%    f is a structure representing a convex function. Inside the structure, there
%     have to be the prox of the function that can be called by f.prox and 
%     the  function itself that can be called by f.eval. 
%
%    A is the operator
%
%    At is the adjoint operator of A
%
%    param a Matlab structure containing solver paremeters. See the
%     function SOLVEP for more information. Additionally it contains those
%     aditional fields:  
%
%      param.nu : bound on the norm of the operator A (default: 1), i.e.
%
%        ` ||A x||^2 <= nu * ||x||^2 
%
%      param.method : is the method used to solve the problem. It can be 'FISTA' or
%       'ISTA'. By default, it's 'FISTA'.
%            
%   See also:  forward_backward solvep admm
%
%   References:
%     P. Combettes and J. Pesquet. A douglas--rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564--574, 2007.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/rlr.php

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
% Date: sept 30 2011
%

% Start the time counter
t1 = tic;

% test the evaluate function
[f] = test_eval(f);

% Optional input arguments
if nargin<5, param=struct; end

if ~isfield(param, 'nu'), param.nu=1 ; end
if ~isfield(param, 'method'), param.method='FISTA' ; end


% setting the function f2 
f2.grad = @(x) 2*At((A(x)-x_0));
f2.eval = @(x) (norm(A(x)-x_0,'fro'))^2;
f2.beta = 2 * param.nu;

[sol, info] = solvep(x_0,{f,f2},param);

info.algo=mfilename;
info.time=toc(t1);


end

