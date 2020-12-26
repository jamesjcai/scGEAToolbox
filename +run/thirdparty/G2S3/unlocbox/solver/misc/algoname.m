function algo = algoname(name)
%ALGONAME return an algorithm from is name
%   Usage: algo = algoname(name)
%
%   Input parameters:
%       name    : name of the algorithm (string)
%   Output parameters:
%       algo    : algorithm (struct)
%
%   The structure algo contains 3 fields:
%    algo.name : the name of the algorithm (string)
%    algo.ignitialize : the initialization function of the algorithm
%    algo.algorithm : the core of one iteration of the algorithm
%    algo.finalize : post processing
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/misc/algoname.php

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


switch lower(name)
    case 'forward_backward'
        algo = forward_backward_alg();
    case 'douglas_rachford'
        algo = douglas_rachford_alg();
    case 'admm'
        algo = admm_alg();
    case 'sdmm'
        algo = sdmm_alg();
    case 'ppxa'
        algo = ppxa_alg();        
    case 'generalized_forward_backward'
        algo = generalized_forward_backward_alg();      
    case 'gradient_descent'
        algo = gradient_descent_alg();
%     case 'backward_backward'
%         algo = backward_backward_alg();  
    case 'pocs'
        algo = pocs_alg();
    case 'chambolle_pock'
        algo = chambolle_pock_alg();      
    case 'fb_based_primal_dual'
        algo = fb_based_primal_dual_alg();      
    case 'fbf_primal_dual'
        algo = fbf_primal_dual_alg();  
    otherwise
        error('Unknown algorithm name')
        
end

end


