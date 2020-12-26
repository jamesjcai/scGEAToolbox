function [ ] = setfig( param )
%SETFIG Set default parameters for plotting
%   Usage: setfig(param);
%   
%   Input parameters:
%       param   : optional parameters
%   Output parameters:
%       none
%
%   param a Matlab structure containing the following fields:
%
%    param.position : position and size of the figure 
%     (default [100 100 600 400])
%    param.labelsize : Size of the label (default 12)
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/misc/setfig.php

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




% Nathanael Perraudin
% 25 November 2013

if nargin<1
    param=struct;
end

% Optional parameters
if ~isfield(param, 'position'), param.position = [100 100 250 250]; end
if ~isfield(param, 'labelsize'), param.labelsize = 12; end



% set the axes
   
%    set(0,'DefaultFigurePaperPosition',param.position)
    set(0,'DefaultFigurePosition',param.position);
    %set(0,'DefaultFigurePaperPosition','auto');
   
    % Change default axes fonts.
    set(0,'DefaultAxesFontSize',  param.labelsize)

    % Change default text fonts.
    set(0,'DefaultTextFontSize', param.labelsize)
    


end


