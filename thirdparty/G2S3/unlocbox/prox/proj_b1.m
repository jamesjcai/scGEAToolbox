function [sol,info] = proj_b1(x, ~, param)
%PROJ_B1 Projection onto a L1-ball
%   Usage:  sol=proj_b1(x, ~, param)
%           [sol,infos]=proj_b1(x, ~, param)
%
%   Input parameters:
%         x     : Input signal.
%         param : Structure of parameters.
%   Output parameters:
%         sol   : Solution.
%         info  : Structure summarizing informations at convergence
%
%   PROJ_B1(x,~,param) solves:
%
%      sol = argmin_{z} ||x - z||_2^2   s.t.  ||w.*z||_1 < epsilon
%
%   Remark: the projection is the proximal operator of the indicative function of
%   w.*z||_1 < epsilon. So it can be written:
%
%      prox_{f, gamma }(x)      where       f= i_c(||w.*z||_1 < epsilon)
%
%   param is a Matlab structure containing the following fields:
%
%    param.epsilon : Radius of the L1 ball (default = 1).
%
%    param.weight : contain the weights (default ones).
%
%    param.verbose : 0 no log, 1 a summary at convergence, 2 print main
%     steps (default: 1)
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
%    info.final_eval : Final evaluation of the function
%
%    info.crit : Stopping critterion used 
%
%
%   Rem: The input "~" is useless but needed for compatibility issue.
%
%   This code is partly borrowed from the SPGL toolbox!
%
%   See also:  proj_b2 prox_l1
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/proj_b1.php

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

%
% Author: Nathanael Perraudin
% Date: February 2015
% Testing: test_proj_b1

% Start the time counter
t1 = tic;

% Optional input arguments
if nargin<3, param=struct; end
if ~isfield(param, 'epsilon'), param.epsilon = 1; end
if ~isfield(param, 'verbose'), param.verbose = 1; end
if ~isfield(param, 'weight'), param.weight = ones(size(x)); end

if isfield(param,'w')
    error('Change in the UNLocBoX! Use weight instead of w!');
end

if isscalar(param.weight), param.weight = ones(size(x))* param.weight; end
param.weight = abs(param.weight);

% Quick return for the easy cases.
if sum(param.weight) == 0
   sol   = x;
   iter = 0;

    crit='--';
    info.algo=mfilename;
    info.iter=iter;
    info.final_eval=0;
    info.crit=crit;
    info.time=toc(t1);
   return
end

% Get sign of b and set to absolute values
signx = sign(x);
x = abs(x);

idx = find(x > eps); % Get index of all non-zero entries of d
sol   = x;             % Ensure x_i = b_i for all i not in index set idx
[sol(idx),iter] = one_projector(sol(idx),param.weight(idx),param.epsilon);


% Restore signs in x
sol = sol.*signx;



% Log after the projection onto the L2-ball
if param.verbose >= 1
    fprintf('  Proj. B1: epsilon = %e, ||x||_2 = %e,\n', param.epsilon, norm(sol,1));
end

crit='--';
info.algo=mfilename;
info.iter=iter;
info.final_eval=norm(param.weight.*sol,1);
info.crit=crit;
info.time=toc(t1);

end




function [sol,iter] = one_projector(x,weight,tau)
% This code is partly borrowed from the SPGL toolbox
  % Initialization
   N = length(x);
   sol = zeros(N,1);

   % Check for quick exit.
   if (tau >= norm(weight.*x,1)), sol = x; iter = 0; return; end
   if (tau <  eps         ),        iter = 0; return; end

   % Preprocessing (b is assumed to be >= 0)
   [sw,idx] = sort(x ./ weight,'descend'); % Descending.
   x  = x(idx);
   weight  = weight(idx);

   % Optimize
   csdb = 0; csd2 = 0;
   soft = 0; ii = 1;
   while (ii <= N)
      csdb = csdb + weight(ii).*x(ii);
      csd2 = csd2 + weight(ii).*weight(ii);
  
      alpha1 = (csdb - tau) / csd2;
      alpha2 = sw(ii);

      if alpha1 >= alpha2
         break;
      end
    
      soft = alpha1;  ii = ii + 1;
   end
   sol(idx(1:ii-1)) = x(1:ii-1) - weight(1:ii-1) * max(0,soft);

   % Set number of iterations
   iter = ii;

end

