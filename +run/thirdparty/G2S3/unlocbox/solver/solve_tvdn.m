function [sol, info] = solve_tvdn(y, epsilon, A, At, param)
%SOLVE_TVDN Solve TVDN problem
%   Usage: sol = solve_tvdn(y, epsilon, A, At, param)
%          sol = solve_tvdn(y, epsilon, A, At)
%          [sol, info] = solve_tvdn(...)
%
%   Input parameters:
%         y     : Measurements
%         epsilon: Radius of the L2 ball
%         A     : Operator
%         At    : Adjoint of A
%         param : Optional parameter
%   Output parameters:
%         sol   : Solution
%         info  : Structure summarizing informations at convergence
%
%   sol = SOLVE_TVDN(Y, epsilon, A, At, PARAM) solves:
%
%      sol arg min ||x||_TV  s.t.  ||y-A x||_2 < epsilon
%
%   Y contains the measurements. A is the forward measurement operator and
%   At the associated adjoint operator. PARAM a Matlab structure containing
%   the following fields:
%
%   General parameters:
% 
%    param.verbose : 0 no log, 1 print main steps, 2 print all steps.
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.useGPU : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
%
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     where  n(t) = (x)||_{TV} is the objective function at iteration t*
%     by default, tol=10e-4.
%
%    param.gamma : control the converge speed (default: 1e-1).
% 
% 
%   Projection onto the L2-ball :
%
%    param.tight_b2 : 1 if A is a tight frame or 0 if not (default = 1)
% 
%    param.nu_b2 : bound on the norm of the operator A, i.e.
%
%        ` ||A x||^2 <= nu * ||x||^2 
%
%    param.tol_b2 : tolerance for the projection onto the L2 ball (default: 1e-3):
%
%        epsilon/(1-tol) <= ||y - A z||_2 <= epsilon/(1+tol)
%
%    param.maxit_b2 : max. nb. of iterations for the projection onto the L2
%     ball (default 200).
% 
% 
%   Proximal TV operator:
%
%    param.maxit_tv : Used as stopping criterion for the proximal TV
%     operator. Maximum number of iterations.
%
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%
%    info.iter : Number of iteration
%
%    info.time : Time of exectution of the function in sec.
%
%    info.final_eval : Final evaluation of the objectivs functions
%
%    info.crit : Stopping critterion used 
%
%    info.rel_norm : Relative norm at convergence 
%
%    info.residue : Final residue 
%
%
%   The problem is solved thanks to a Douglas-Rachford splitting
%   algorithm.
%
%   Demos: demo_tvdn
%
%   References:
%     P. Combettes and J. Pesquet. A douglas--rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564--574, 2007.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/solve_tvdn.php

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

% Author: Gilles Puy, Nathanael Perraudin
% Date: Nov. 1, 2012


% Optional input arguments
if nargin<5, param=struct; end


% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'tol'), param.tol = 1e-4; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
%if ~isfield(param, 'gamma'), param.gamma = 1e-2; end
if ~isfield(param, 'useGPU'), param.useGPU = 0; end

% Input arguments for projection onto the L2 ball
param_b2.A = A; param_b2.At = At;
param_b2.y = y; param_b2.epsilon = epsilon;
param_b2.verbose = param.verbose-1;
if isfield(param, 'nu_b2'), param_b2.nu = param.nu_b2; end
if isfield(param, 'tol_b2'), param_b2.tol = param.tol_b2; end
if isfield(param, 'tight_b2'), param_b2.tight = param.tight_b2; end
if isfield(param, 'maxit_b2')
    param_b2.maxit = param.maxit_b2;
end

% Input arguments for prox TV
param_tv.verbose = param.verbose-1; param_tv.tol = param.tol;
param_tv.useGPU = param.useGPU;
if isfield(param, 'maxit_tv')
    param_tv.maxit= param.maxit_tv;
end

f1.prox = @(x,T) prox_tv(x,T,param_tv);
f1.eval = @(x) norm_tv(x);

f2.prox = @(x,T) proj_b2(x,T,param_b2);
f2.eval = @(x) eps;

[sol,info] = douglas_rachford(At(y), f2, f1, param);


end

