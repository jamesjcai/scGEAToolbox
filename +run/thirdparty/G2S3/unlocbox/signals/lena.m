function [im] = lena(color)
%LENA  Load the 'lena' test signal
%   Usage: im = lena();
%          im = lena(color);
%
%   Input parameters:
%       color   : boolean 
%   Output parameters:
%       none    :
%
%   LENA() loads the graylevel 'lena' signal. Lena is a common
%   image processing test image of resolution (512 x 512). However, we do
%   recommand not to use it. 
%
%   LENA(1) loads the color 'lena' signal.
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
%       im = lena();
%       imagescgray(im);
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/signals/lena.php

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
  
if nargin < 1
  color = 0;
end;

f = mfilename('fullpath');

% Load the signal

im = imread([f, '.png']);

if ~ color
    im = rgb2gray(im);
end

im = double(im) / 255;

