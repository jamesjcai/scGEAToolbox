function snr_val = snr(map_init, map_noisy)
%SNR Compute the SNR between two maps
%   Usage snr_val = snr(map_init, map_noisy) 
% 
%   Input parameters:
%         map_init : initial signal
%         map_recon: noisy signal
%   Output parameters:
%         snr_val  : snr
%
%   computes the SNR between the maps map_init
%   and map_noisy. The SNR is computed as:
%
%       10  log10( var(map_init) / var(map_init-map_noisy) )
%
%   where var stands for the matlab built-in function that computes the
%   variance.
% 
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/snr.php

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
 
% Author: Gilles Puy
% Date: 2009
% 

noise = map_init(:)-map_noisy(:);
var_init = var(map_init(:));
var_den = var(noise(:));
snr_val = 10 * log10(var_init/var_den);

end

