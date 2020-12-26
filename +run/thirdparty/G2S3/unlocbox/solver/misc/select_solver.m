function solver = select_solver(Fgrad, Fprox)
%SELECT_SOLVER This function choose a default solver
%   Usage: solver = select_solver(Fg, Fp)
%
%   This function chooses a default solver depending on the function Fg and
%   Fp. Fg is cell array of smooth functions and Fp a cell array of non
%   smooth functions.
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/misc/select_solver.php

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

    n_linear = numberofL(Fprox);
    n_prox = numel(Fprox);
    n_grad = numel(Fgrad);
    
    if n_grad > 0       % There are smooth functions
        if n_linear > 2
            error('Sorry, no solver is able to solve your problem yet!')
        end
        if n_prox == 0
            solver = 'GRADIENT_DESCENT';
        elseif n_prox == 2 
            solver = 'FB_BASED_PRIMAL_DUAL';
        elseif (n_prox <= 2) && (n_linear > 0)
            solver = 'FB_BASED_PRIMAL_DUAL';            
        elseif n_prox == 1
            solver = 'FORWARD_BACKWARD';
        else
            solver = 'GENERALIZED_FORWARD_BACKWARD';
        end
    else    % All function are non smooth
     
        if (n_prox <= 2) && (n_linear == 1)
            solver = 'CHAMBOLLE_POCK';    
        elseif (n_linear > 1)
            solver = 'SDMM';
        elseif n_prox == 2
            solver = 'DOUGLAS_RACHFORD';
        else
            solver = 'PPXA';
        end
    end
        
end




function n = numberofL(Fp)
% Return the number of functions with a linear operator inside
n = 0;
    for ii = 1:length(Fp)
        if isfield(Fp{ii},'L') 
            n = n + 1;
            if ~isfield(Fp{ii},'Lt')
                warning('You did not define the Lt operator!');
            end
        end
    end
end
