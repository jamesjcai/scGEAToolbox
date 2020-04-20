function gamma = log_decreasing_ts(x, gamma_in, gamma_fin, nit)
%LOG_DECREASING_TS Log decreasing timestep for UNLCOBOX algorithm
%   Usage gamma = log_decreasing_ts(x, gamma_in, gamma_fin, nit);
%
%   Input parameters:
%         x         : Structure of data
%         gamma_in  : Initial timestep
%         gamma_fin : Final timestep
%         nit       : Number of iteration for the decrease
%
%   Output parameters:
%         gamma     : Timestep at iteration t
%
%   This plug-in computes a new timestep at each iteration. It makes a log
%   decreasing timestep from gamma_in to gamma_fin in nit iterations.
%   To use this plugin, define:
%
%       param.do_ts = @(x) log_decreasing_ts(x, gamma_in, gamma_fin, nit);
%
%   in the structure of optional argument of the solver.
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/plugins/log_decreasing_ts.php

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
% Date:   3 April 2014


    if x.iter > nit
        gamma = gamma_fin;
    else
        ts = gamma_in./linspace(1,gamma_in/gamma_fin,nit);
        gamma = ts(x.iter);
    end
end

