function [ sol ] = plot_objective(info_iter, fig)
%PLOT_OBJECTIVE Plot objective function over iters(plugin for the UNLocBoX)
%   Usage [ sol ] = plot_objective( info_iter, fig );
%
%   Input parameters:
%         info_iter   : Structure of info of current iter of algorithm
%         fig   : Figure
%
%   Output parameters:
%         sol   : Current solution
%
%   This plugin displays the image every iterations of an algorithm. To use
%   the plugin juste define:
%       
%       fig = figure(100);
%       param.do_sol = @(x) plot_objective(x, fig);
%
%   In the structure of optional argument of the solver.
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/plugins/plot_objective.php

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
% Date  : November 2013




% select the figure
if info_iter.iter<2
    figure(fig);
end

% 
title(['Current it: ', num2str(info_iter.iter),'   Curr obj: ', ...
    num2str(info_iter.info.objective(end))]);
semilogy(info_iter.info.objective); title('Objective function')
drawnow;

% return the image
sol=info_iter.sol;

end


