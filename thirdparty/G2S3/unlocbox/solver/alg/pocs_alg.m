function s = pocs_alg()
   s.name = 'POCS';
   s.initialize = @(x_0, fg, Fp, param) pocs_initialize(x_0,fg,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) pocs_algorithm(Fp, sol, s);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = pocs_initialize(x_0,fg,param)
    


    sol = x_0;
    s = {};
    if fg.beta
        error('Beta = 0! This solver requires only projections functions.');
    end
end


function [sol, s] = pocs_algorithm(Fp, sol, s)
    
    for ii = 1 : length(Fp)
       sol = Fp{ii}.prox(sol,0);
    end

end

%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/alg/pocs_alg.php

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

