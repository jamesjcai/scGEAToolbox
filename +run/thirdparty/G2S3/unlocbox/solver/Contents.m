% Unlocbox - Solvers
%
%  Universal solver
%    SOLVEP             -  Universal solver
%
%  General solvers
%    ADMM               -  Alternating-direction method of multipliers algorithm.
%    DOUGLAS_RACHFORD   -  Foward backward splitting algorithm.
%    FORWARD_BACKWARD   -  Foward backward splitting algorithm.
%    GENERALIZED_FORWARD_BACKWARD - Generaliyed foward backward splitting algorithm.
%    PPXA               -  Parallel Proximal algorithm.
%    SDMM               -  Simultaneous-direction method of multipliers algorithm.
%    FB_BASED_PRIMAL_DUAL -  Forward Backard based Primal Dual
%    FBF_PRIMAL_DUAL    -  Forward backward forward primal dual
%    GRADIENT_DESCENT   -  Simple gradient descent solver.
%    POCS               -  Projection onto convex sets.
%
%  Composed solvers
%    RLR                -  Regularized Linear Regresssion solver (special case of forward-backward).
%    SOLVE_BPDN         -  Solve a BPDN (Basis Pursuit denoising) problem.
%    SOLVE_TVDN         -  Solve a TVDN (TV denoising) problem.
%
%  Demo solver
%    DEMO_FORWARD_BACKWARD_ALG - Example solver
%
%  For help, bug reports, suggestions etc. please send an email to
%  unlocbox (at) groupes (dot) epfl (dot) ch
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/Contents.php

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
% see also: solvep prox
