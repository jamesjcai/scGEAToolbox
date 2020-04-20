function [im]=barbara()
%BARBARA  Load the 'barbara' test signal
%
%   BARBARA loads the 'barbara' signal. Barbara is an image commonly
%   used in image compression and filtering papers because it contains a
%   range of tones and many thin line patterns. The resolution is (512 x
%   512). 
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
%       im = barbara();
%       imagescgray(im);
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/signals/barbara.php

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
  
if nargin>0
  error('This function does not take input arguments.')
end;

f=mfilename('fullpath');

% Load the signal
im = double(imread([f,'.png']))/255;


