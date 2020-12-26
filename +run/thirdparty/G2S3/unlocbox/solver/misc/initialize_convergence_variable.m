function [info, iter, s] = initialize_convergence_variable(sol, s, fg, Fp, param)
%INITIALIZE_CONVERGENCE_VARIABLE Initialize all convergence variable
%
%   Handle the initialization of info depending on the stopping criterion.
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/misc/initialize_convergence_variable.php

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
% Date: 16 oct 2015

% Test if the the dual variables are present
if strcmpi(param.stopping_criterion, 'rel_norm_primal_dual') ...
    || strcmpi(param.stopping_criterion, 'rel_norm_dual');
    if ~isfield(s,'dual_var')
        error('The dual variable are not ignitialized. Put them in s.dual_var');
    end
end


% In debut mode we inialize everything
if param.debug_mode
    info.objective = nan(param.maxit+1, 1);
    info.objective(1) = eval_function(fg, Fp, sol, s, param);
    s.prev_sol = sol;
    info.rel_norm_primal = nan(param.maxit, 1);
    info.rel_norm_dual = nan(param.maxit, 1);
    info.rel_eval = nan(param.maxit,1);
    if isfield(s,'dual_var')
    	s.dual_var_old = s.dual_var;
    end
else

    switch lower(param.stopping_criterion)
        case 'rel_norm_obj'
            info.objective = nan(param.maxit+1,1);
            info.rel_eval = nan(param.maxit,1);
            info.objective(1) = eval_function(fg, Fp, sol, s, param);
        case 'rel_norm_primal'
            s.prev_sol = sol;
            info.rel_norm_primal = nan(param.maxit, 1);
        case 'rel_norm_dual'
            s.dual_var_old = s.dual_var;
            info.rel_norm_dual = nan(param.maxit, 1);
        case 'rel_norm_primal_dual'
            s.dual_var_old = s.dual_var;
            s.prev_sol = sol;
            info.rel_norm_primal = nan(param.maxit, 1);
            info.rel_norm_dual = nan(param.maxit, 1);           
        case 'obj_increase'
            info.objective = nan(param.maxit+1,1);
            info.objective(1) = eval_function(fg, Fp, sol, s, param);
        case 'obj_threshold'
            info.objective = nan(param.maxit+1,1);
            info.objective(1) = eval_function(fg, Fp, sol, s, param);   
        otherwise
            error('Unknown stopping criterion!')
    end
end   


% We start at iteration 1
iter = 1;


end
