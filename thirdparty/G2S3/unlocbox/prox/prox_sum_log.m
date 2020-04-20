function [sol,info] = prox_sum_log(x, gamma, param)
%PROX_SUM_LOG Proximal operator of log-barrier  - sum(log(x))
%   Usage:  sol = prox_sum_log(x, gamma)
%           sol = prox_sum_log(x, gamma, param)
%           [sol, info] = prox_sum_log(x, gamma, param)
%
%   Input parameters:
%         x     : Input signal (vector or matrix!).
%         gamma : Regularization parameter.
%         param : Structure of optional parameters.
%
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   PROX_SUM_LOG(x, gamma, param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 - gamma * sum(log(z))
%
%   param is a Matlab structure containing the following fields:
%
%    param.verbose : 0 no log, (1) print -sum(log(z)), 2 additionally
%     report negative inputs.
%   
%   MATRICES:
%   Note that this prox works for matrices as well. The log of the sum
%   gives the same result independently of which dinension we perform the
%   summation over:
%
%       sol = (x + sqrt(x.^2 + 4*gamma)) /2;
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%    info.iter : Number of iteration
%    info.time : Time of exectution of the function in sec.
%    info.final_eval : Final evaluation of the function
%    info.crit : Stopping critterion used 
%
%
%   See also:  prox_l1 prox_l2 prox_tv prox_sum_log_norm2
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/prox_sum_log.php

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
% Date: June 2015
% Testing: test_prox_sum_log

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



sol = (x + sqrt(x.^2 + 4*gamma)) /2;
info.algo = mfilename;
info.iter = 0;
info.final_eval = -gamma * sum(log(x(:)));
info.crit = '--';
info.time = toc(t1);
    
% Log after the prox
if param.verbose >= 1
    fprintf('  prox_sum_log: - sum(log(x)) = %e', info.final_eval / gamma);
    if param.verbose > 1
        n_neg = nnz(x(:) <= 0);
        if n_neg > 0
            fprintf(' (%d negative elements, log not defined, check stability)', n_neg);
        end
    end
    fprintf('\n');
end


end



