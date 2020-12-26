function [sol,info] = prox_log_det(x, gamma, param)
%PROX_LOG_DET Proximal operator of log-determinant
%   Usage:  sol = prox_log_det(x, gamma)
%           sol = prox_log_det(x, gamma, param)
%           [sol, info] = prox_log_det(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal.
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   prox_l1(x, gamma, param) solves:
%
%      sol = argmin_{Z} 0.5*||X - Z||_2^2 - gamma * log(det(Z))
%
%   param is a Matlab structure containing the following fields:
%
%    param.verbose : 0 no log, (1) print -(log(det(z)), 2 additionally
%       report negative inputs.
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
%   See also:  prox_l1 prox_l2 prox_tv
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/prox_log_det.php

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

% Author: Vassilis Kalofolias
% Date: September 2015
% Testing: TODO!!  test_prox_log_det 

% Start the time counter
t1 = tic;

if nargin < 3, param = struct; end

if ~isfield(param, 'verbose'), param.verbose = 1; end


% test the parameters
% stop_error = test_gamma(gamma);
if gamma < 0
    error('Gamma can not be negative');
elseif gamma == 0
    stop_error = 1;
else
    stop_error = 0;
end

if stop_error
    sol = x;
    info.algo = mfilename;
    info.iter = 0;
    info.final_eval = 0;
    info.crit = '--';
    info.time = toc(t1);
    return
end


[Q, e] = eig(x, 'vector');
sol = Q*diag(e + sqrt(e.^2 + 4*gamma)/2)*Q';
info.algo = mfilename;
info.iter = 0;
info.final_eval = -gamma * sum(log(e));
info.crit = '--';
info.time = toc(t1);
    
% Log after the prox
if param.verbose >= 1
    fprintf('  prox_log_det: - log(det(x)) = %e', info.final_eval / gamma);
    if param.verbose > 1
        n_neg = nnz(e(:) <= 0);
        if n_neg > 0
            fprintf(' (%d negative eigenvalues, log det not defined, check stability)', n_neg);
        end
    end
    fprintf('\n');
end


end



