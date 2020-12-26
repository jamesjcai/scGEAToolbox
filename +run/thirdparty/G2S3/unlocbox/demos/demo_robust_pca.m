%DEMO_ROBUST_PCA
%
%      argmin_Z || Z - Zn ||_1 + tau || Z ||_*
%
%   where tau is the regularization parameter.
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/demos/demo_robust_pca.php

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

%% parameter
verbose = 2; % verbosity parameter

N = 100; % number of columns
M = 200; % number of rows
k = 10; % rank
p = 0.1; % probability for the sparse noise 
sigma = 10; % noise level
tau = 10; % regularization parameter
%% Create the data
X = randn(100,k);
Y = rand(k,M);


Z = X*Y; % Low rank data
N = sigma*sprand(N,M,p);
Znoisy = Z + N; % measurements

%% Define the function inside the problem


paraml1.verbose = verbose-1;
paraml1.y = Znoisy;
f_f.prox = @(x,T) prox_l1(x,T,paraml1);
f_f.eval = @(x) sum(abs(x(:)));

% 2) nuclear norm f_n (X) = tau || X ||_*
paramnuclear.verbose = verbose -1;
f_n.prox = @(x,T) prox_nuclearnorm(x,T*tau,paramnuclear);
f_n.eval = @(x) tau*norm_nuclear(x);



%% Solve the problem
paramsolver.verbose = verbose;
paramsolver.gamma = 0.5; % timestep.

Zsol = solvep(Znoisy, {f_f,f_n}, paramsolver);

N = Znoisy - Zsol;

%% Errors + plots

err_sol = norm(Zsol-Z,'fro')/norm(Z,'fro')
err_in = norm(Znoisy-Z,'fro')/norm(Z,'fro')

figure(1)
subplot(221)
imagesc(Z)
title('Original low rank data')
subplot(222)
imagesc(Znoisy)
title('Noisy data')
subplot(223)
imagesc(Zsol)
title('Recovery')

subplot(224)
imagesc(N)
title('Sparse noise')
