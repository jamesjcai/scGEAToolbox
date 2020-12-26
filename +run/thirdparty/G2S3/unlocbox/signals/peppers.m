function [im]=peppers(color)
%PEPPERS  Load the 'peppers' test signal
%   Usage: im = peppers();
%          im = peppers(color);
%
%   Input parameters:
%       color   : boolean 
%   Output parameters:
%       none    :
%
%   PEPPERS() loads the graylevel 'peppers' signal. Peppers is a common
%   image processing test image of resolution (512 x 512).   
%
%   PEPPERS(1) loads the color 'peppers' signal.
%   
%   This signal, and other standard image tests signals, can be found on
%   Morgan McGuire's Computer Graphics
%   Archivehttp://graphics.cs.williams.edu/data/images.xml. 
%
%   For convenience the output image is normalized by 255 and converted to
%   double.
%
%   Example
%   -------
%   
%   Load the image and display it:
%
%       im = peppers();
%       imagescgray(im);
%  
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/signals/peppers.php

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
% Date: 25 November 2013
  
if nargin<1
  color = 0;
end;

f=mfilename('fullpath');

% Load the signal

im = imread([f,'.png']);

if ~color
    im = rgb2gray(im);
end

im = double(im)/255;


