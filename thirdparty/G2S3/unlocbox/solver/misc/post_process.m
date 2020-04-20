function [sol, param] = post_process(sol, iter, info, param)
%POST_PROCESS Post processing for the UNLocBoX
%   Usage: [sol,param]=post_process(sol,iter,curr_norm,prev_norm,param);
%
%   
%   This function make evaluate the post processing job in the UnLocBoX
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/misc/post_process.php

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

if ~isfield(param, 'gamma'), param.gamma = 1; end;
if ~isfield(param, 'do_sol'), param.do_sol = @(x) x.sol; end
if ~isfield(param, 'do_ts'), param.do_ts = @(x) x.gamma ; end

info_iter.sol = sol;
info_iter.iter = iter;
info_iter.gamma = param.gamma;
info_iter.info = info;

% Perform tunable user actions every iterations like display the current
% solution or update the timesteps
sol = param.do_sol(info_iter);
param.gamma = param.do_ts(info_iter);



end


