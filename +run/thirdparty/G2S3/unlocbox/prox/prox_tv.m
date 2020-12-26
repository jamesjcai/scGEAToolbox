function [sol, info] = prox_tv(b, gamma, param)
%PROX_TV Total variation proximal operator
%   Usage:  sol=prox_tv(x, gamma)
%           sol=prox_tv(x, gamma,param)
%           [sol, info]=prox_tv(...)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   This function compute the 2 dimentional TV proximal operator evaluated
%   in b. If b is a cube, this function will evaluate the TV proximal
%   operator on each image of the cube. For 3 dimention TV proximal
%   operator the function prox_tv3d can be used.
%
%   PROX_TV(y, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 + gamma * ||z||_TV
%
%   param is a Matlab structure containing the following fields:
%   
%    param.tol : is stop criterion for the loop. The algorithm stops if
%
%         (  n(t) - n(t-1) )  / n(t) < tol,
%      
%     where  n(t) = f(x)+ 0.5 X-Z_2^2 is the objective function at iteration t*
%     by default, tol=10e-4.
%
%    param.maxit : max. nb. of iterations (default: 200).
%
%    param.useGPU : Use GPU to compute the TV prox operator. Please prior 
%     call init_gpu and free_gpu to launch and release the GPU library (default: 0).
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%    param.weights : weights for each dimention (default [1, 1])
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%
%    info.iter : Number of iteration
%
%    info.time : Time of exectution of the function in sec.
%
%    info.final_eval : Final evaluation of the function
%
%    info.crit : Stopping critterion used 
%
%
%   See also:  prox_l1 prox_tv3d prox_tv1d gradient_op div_op
%
%
%   References:
%     A. Beck and M. Teboulle. Fast gradient-based algorithms for constrained
%     total variation image denoising and deblurring problems. Image
%     Processing, IEEE Transactions on, 18(11):2419--2434, 2009.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/prox_tv.php

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


% Author: Nathanael Perraudin, Gilles Puy, Eyal Hirsch
% Date: Jan 2013
%

% Start the time counter
t1 = tic;

% for the GPU
global GLOBAL_useGPU; 

if ~size(GLOBAL_useGPU,1), GLOBAL_useGPU = 0; end

% Optional input arguments

if nargin<3, param=struct; end

if ~isfield(param, 'tol'), param.tol = 10e-4; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'maxit'), param.maxit = 200; end
if ~isfield(param, 'weights'), param.weights = [1, 1]; end
if ~isfield(param, 'useGPU')
    param.useGPU = (GLOBAL_useGPU) || (isa(b,'gpuArray')); 
end

% Test of gamma
if test_gamma(gamma)
    sol = b;
    info.algo=mfilename;
    info.iter=0;
    info.final_eval=0;
    info.crit='--';
    info.time=toc(t1);
    return; 
end

if param.useGPU
    %gpuDevice(1);
    gamma=gpuArray(gamma);
    if isa(b,'gpuArray')
        allGPU=1;
    else
        b=gpuArray(b);
        allGPU=0;
    end
    % Initializations
    [r, s] = gradient_op(b*0);
    pold = r; qold = s;
    told = gpuArray(1); prev_obj = gpuArray(0); 
    verbose=gpuArray(param.verbose);
    tol=gpuArray(param.tol);
else
    % Initializations
    [r, s] = gradient_op(b*0);
    pold = r; qold = s;
    told = 1; prev_obj = 0;
    verbose=param.verbose;
    tol=param.tol;
end


wx = param.weights(1);
wy = param.weights(2);
mt = max(param.weights);

% Main iterations
if verbose > 1
    fprintf('  Proximal TV operator:\n');
end


    
    
for iter = 1:param.maxit

    % Current solution
    sol = b - gamma*div_op(r, s, wx, wy);

    % Objective function value
    tmp = gamma * sum(norm_tv(sol, wx, wy));
    obj = .5*norm(b(:)-sol(:), 2)^2 + tmp;
    rel_obj = abs(obj-prev_obj)/obj;
    prev_obj = obj;

    % Stopping criterion
    if verbose>1
        fprintf('   Iter %i, obj = %e, rel_obj = %e\n', ...
            iter, obj, rel_obj);
    end
    if rel_obj < tol
        crit = 'TOL_EPS'; break;
    end

    % Udpate divergence vectors and project
    [dx, dy] = gradient_op(sol, wx, wy);
    
    r = r - 1/(8*gamma)/mt^2 * dx; 
    s = s - 1/(8*gamma)/mt^2 * dy;
    
    weights = max(1, sqrt(abs(r).^2+abs(s).^2));
    
    p = r./weights; 
    q = s./weights;

    % FISTA update
    t = (1+sqrt(4*told.^2))/2;
    r = p + (told-1)/t * (p - pold); pold = p;
    s = q + (told-1)/t * (q - qold); qold = q;
    told = t;

end


% Log after the minimization
if ~exist('crit', 'var'), crit = 'MAX_IT'; end

if verbose >= 1
    if param.useGPU
        fprintf(['  GPU Prox_TV: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    else
        fprintf(['  Prox_TV: obj = %e, rel_obj = %e,' ...
            ' %s, iter = %i\n'], obj, rel_obj, crit, iter);
    end
end



if param.useGPU
    if ~allGPU
        sol=gather(sol);
    end
    info.iter=gather(iter);
    info.final_eval=gather(obj);
else
    info.iter=iter;
    info.final_eval=obj;
end

info.algo=mfilename;
info.crit=crit;
info.final_eval = tmp;
info.time=toc(t1);

end

