function [sol, info] = proj_box(x, ~, param)
%PROJ_BOX Projection onto the box set (multidimensional interval constraint)
%
%   Usage:  sol = proj_box(x, [])
%           sol = proj_box(x)
%           sol = proj_box(x, [], param)
%           [sol, info] = proj_box(x, [], param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of optional parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing information at convergence
%
%   prox_box(x, [], param) solves:
%
%      sol = argmin_{z} 0.5*||x - z||_2^2 s.t. z < zmax and z > zmin
%
%   where zmax and zmin might be scalar or vector valued.
%
%   param is a Matlab structure containing the following fields:
%
%    param.lower_lim : lower bound(s) for z (default 0)
%    param.upper_lim : upper bound(s) for z (default 1)
%
%   if these two are vector-valued, bounds apply entry-by-entry
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
%
%
%   info is a Matlab structure containing the following fields:
%
%    info.algo : Algorithm used
%
%    info.iter : Number of iterations (this function is not iterative)
%
%    info.time : Time of exectution of the function in sec.
%
%    info.final_eval : Final evaluation of the function
%
%    info.crit : Stopping critterion used (one shot here)
%
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   See also:  proj_b2
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/proj_box.php

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


% Author: Pavel Rajmic
% Date: June 2015
%




% Start the time counter
t1 = tic;

%% Optional input arguments
if nargin<3, param=struct; end

% Optional input arguments
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'lower_lim'), param.lower_lim = 0; end
if ~isfield(param, 'upper_lim'), param.upper_lim = 1; end

%% Check input
if ~isvector(x) | ~isvector(param.lower_lim) | ~isvector(param.upper_lim)
    error('The inputs must be scalars or vectors.')
end

sx = size(x); %size (for future reshape)
lx = max(sx); %length

%make columns
x = x(:);
param.lower_lim = param.lower_lim(:);
param.upper_lim = param.upper_lim(:);

if ~(lx == length(param.lower_lim)) || ~(lx == length(param.upper_lim))
    if ~(length(param.lower_lim) == 1) && ~(length(param.upper_lim) == 1)
        error('Sizes not compatible. See the documentation.')
    end
end

if ~(length(param.lower_lim) == length(param.upper_lim))
    error('Lower and upper constraints must be of same length.')
end

% Check feasibility
if any(param.upper_lim < param.lower_lim)
    error('The feasible set is empty.');
end


%% Projection
sol = x; %copy

greater = sol > param.upper_lim;
lower = sol < param.lower_lim;
    
if length(param.upper_lim) > 1 %vector limits
    %substitute greater values by upper limits
    sol(greater) = param.upper_lim(greater);
    %substitute lower values by lower limits
    sol(lower) = param.lower_lim(lower);
    
else %scalar limit
    %substitute greater values by upper limit
    sol(greater) = param.upper_lim;
    %substitute lower values by lower limit
    sol(lower) = param.lower_lim;
end

sol = reshape(sol, sx);

%% Rest
norm_l2 = 0.5 * norm(x(:)-sol(:))^2;

% Log after the projection
if param.verbose >= 1
    fprintf('  proj_box: 0.5*|| x - z ||_2^2 = %e \n', norm_l2);
end

info.algo = mfilename;
info.iter = 0;
info.final_eval = norm_l2;
info.crit = 'ONE-SHOT';
info.time = toc(t1);

end
