function init_unlocbox()
%INIT_UNLOCBOX Initialize the toolbox
%   Usage: init_unlocbox()
%
%   Initialisation script for the unlocbox
%   This script add the different path needed to run the toolbox
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/init_unlocbox.php

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
% Date: nov 2012


%% adding dependency
global GLOBAL_useGPU;
global GLOBAL_path;
GLOBAL_path = fileparts(mfilename('fullpath'));
GLOBAL_useGPU = 0;

addpath(genpath(GLOBAL_path));

% Load the version number
bp=[GLOBAL_path,filesep];
[FID, MSG] = fopen ([bp,'unlocbox_version'],'r');
if FID == -1
    error(MSG);
else
    unlocbox_version = fgetl (FID);
    fclose(FID);
end

banner = sprintf(strcat(... 
'UnLocBoX version %s. Copyright 2012-2015 LTS2-EPFL, by Nathanael Perraudin'), ...
                   unlocbox_version);
% display banner
disp(banner);

if GLOBAL_useGPU && gpuDeviceCount
    dev=gpuDevice(1);
    if dev.DeviceSupported
        reset(gpuDevice);
        disp('GPU loaded');
    else
        disp(['GPU not loaded.  To remove the previous warning, '...
        'set GLOBAL_useGPU to 0']);  
        GLOBAL_useGPU=0;
    end
end

kbstop('stop');
