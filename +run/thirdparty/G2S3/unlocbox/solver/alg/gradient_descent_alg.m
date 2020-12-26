function s = gradient_descent_alg()
   s.name = 'GRADIENT_DESCENT';
   s.initialize = @(x_0, fg, Fp, param) gradient_descent_initialize(x_0,fg,Fp,param);
   s.algorithm = @(x_0, fg, Fp, sol, s, param) gradient_descent_algorithm(fg, sol, s, param);
   s.finalize = @(x_0, fg, Fp, sol, s, param) sol;

end

function [sol, s, param] = gradient_descent_initialize(x_0,fg,Fp,param)
    

%     s.u_n = x_0;
%     s.tn = 1;
%     sol = s.u_n - param.gamma*fg.grad(s.u_n);
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/alg/gradient_descent_alg.php

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
    sol = x_0;
    s = {};
    
    if numel(Fp)>0
        error('This solver can not be used to optimize non smooth functions')
    end
    
    if ~fg.beta
        error('Beta = 0! This solver requires a smooth term.');
    end
end


function [sol, s] = gradient_descent_algorithm(fg, sol, s, param)
    
    
%     tn1 = (1 + sqrt(1+4*s.tn^2))/2;
%     tmp = sol + (s.tn-1)/tn1*(s.u_n - sol);
%     s.u_n = sol;
%     sol = tmp - param.gamma*fg.grad(tmp);
%     s.tn = tn1;
%     
    sol = sol-param.gamma*fg.grad(sol);


end

