function [ s ] = plot_signal( x,p,fig )
%PLOT_SIGNAL Plot image plugin for the UNLocBoX
%   Usage [ s ] = plot_signal( im,fig );
%
%   Input parameters:
%         x     : Structure of data
%         p     : Number of iteration between 2 plots...
%         fig   : Figure
%
%   Output parameters:
%         s     : Input image
%
%   This plugin display the signal every iterations of an algorithm. To use
%   the plugin juste define:
%       
%       fig=figure(100);
%       param.do_sol=@(x) plot_signal(x,p,fig);
%
%   In the structure of optional argument of the solver.
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/plugins/plot_signal.php

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
% Date  : 3rd april 2014

if ~mod(x.iter-1,p)
    % select the figure
    if x.iter<2
        figure(fig);
    end
    % display the signal
    plot(x.sol);
    title(['Current it: ', num2str(x.iter),'   Curr obj: ', ...
        num2str(x.info.objective(x.iter))]);

    drawnow;
end

% return the signal
s=x.sol;


end


