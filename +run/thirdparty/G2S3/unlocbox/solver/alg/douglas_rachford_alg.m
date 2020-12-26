function s = douglas_rachford_alg()
   s.name = 'DOUGLAS_RACHFORD';
   s.initialize = @(x_0, fg, Fp, param) douglas_rachford_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) douglas_rachford_algorithm(Fp, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = douglas_rachford_initialize(x_0,fg,Fp,param)

    if ~isfield(param, 'lambda'), param.lambda=1 ; end
    s.lambda = param.lambda;
    s.x_n = {};
    s.u_n = x_0;
    sol = x_0;
    if fg.beta
        error('Douglas rachford needs only function with proximal operators')
    end
    
    if ~(numel(Fp)==2)
        error('Douglas rachford needs exactly 2 functions')
    end
    
end


function [sol, s] = douglas_rachford_algorithm(Fp, sol, s, param)
%     s.x_n{1} = Fp{1}.prox_ad(2*sol-s.u_n,param.gamma);
%     s.u_n=s.u_n+param.lambda*(s.x_n{1}{1}-sol);
%     s.x_n{2} =Fp{2}.prox_ad(s.u_n,param.gamma);
%     sol = s.x_n{2}{1};
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/alg/douglas_rachford_alg.php

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
    s.x_n{1} = Fp{1}.prox(2*sol-s.u_n,param.gamma);
    s.u_n=s.u_n+param.lambda*(s.x_n{1}-sol);
    s.x_n{2} =Fp{2}.prox(s.u_n,param.gamma);
    sol = s.x_n{2};
       
end

