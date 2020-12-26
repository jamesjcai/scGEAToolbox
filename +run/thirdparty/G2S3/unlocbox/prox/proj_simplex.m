function [sol, info] = proj_simplex(x, ~, param)
%PROJ_SIMPLEX projection onto the simplex sum(x) = c, x >= 0, c = scalar
%   Usage:  sol = proj_simplex(x, ~, param)
%           [sol, info] = proj_simplex(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info : Structure summarizing informations at convergence
%
%   PROJ_SIMPLEX(x,~,param) solves:
%
%      sol = argmin_{z} ||x - z||_2^2   s.t.  sum(x) = c, x >= 0
%
%   If x is a matrix, the projection is done for each column separately,
%   or across param.dim if specified.
%
%   param is a Matlab structure containing the following fields:
%
%    param.c : scalar (default: 1).
%
%    param.dim : dimension of summation if x is a matrix (default: 1)
%
%    param.verbose : 0 no log, 1 a summary at convergence (default: 0)
%
%   
%   info is a Matlab structure containing the following fields:
%
%    info.time : Time of execution of the function in sec.
%
%    info.final_eval : Final evaluation of the function
%
%
%
%
%   Rem: The input "~" is useless but needed for compatibility reasons.
%
%   See also:  proj_linear_eq, proj_b1
%
%   References:
%     J. Duchi, S. Shalev-Shwartz, Y. Singer, and T. Chandra. Efficient
%     projections onto the l 1-ball for learning in high dimensions. In
%     Proceedings of the 25th international conference on Machine learning,
%     pages 272--279. ACM, 2008.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/proj_simplex.php

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

% TODO: implement "figure 2" algorithm of same paper!! (difficult to
% vectorize!)
%
% Author: Vassilis Kalofolias
% Date: February 2016
% Testing: test_proj_simplex

% Start the time counter
t1 = tic;

if nargin < 3
    param = struct;
end


if ~isfield(param, 'c'), param.c = 1; end
if ~isfield(param, 'dim')
    if isvector(x)
        % 1 for column, 2 for row
        param.dim = 1 + isrow(x);
    else
        param.dim = 1;
    end
end
if ~isfield(param, 'verbose'), param.verbose = 0; end


[m, n] = size(x);

% implementing the algorithm of figure 1 of [1]
x_s = sort(x, param.dim, 'descend');
if param.dim == 1
    one_over_j = 1./(1:m)';
else
    one_over_j = 1./(1:n);
end
    
Thetas = bsxfun(@times, one_over_j, (cumsum(x_s, param.dim)-param.c)  );
% find number of points that are below 0
rho = sum((x_s - Thetas) > 0, param.dim);

if param.dim == 1
    theta = Thetas(sub2ind([m, n], rho, 1:n));
else
    theta = Thetas(sub2ind([m, n], (1:m)', rho));
end
sol = max(0, bsxfun(@plus, x, -theta));

% info about algorithm
info.algo = mfilename;
info.iter = 1;
info.final_eval = 0;
info.crit = '';
info.time = toc(t1);

if param.verbose >= 1
    fprintf('  Proj. simplex:  %e  seconds\n', info.time);
end

end



