%DEMO_COMPRESS_SENSING3 Compress sensing example using grouped L12 norm
%   
%   We present a compress sensing example solved with the douglas rachford
%   solver. The particularity of this example is the use of a mixed norm. We
%   do not only know the the signal is sparse, we also know that the
%   sparse coefficients are grouped.
%
%   The problem can be expressed as this
%
%        argmin || x||_{2,1} s.t ||b-Ax||_2 < epsilon
%
%   Where b are the measurements and A the measurement matrix.
%
%   We set 
%
%    f_1(x)=||x||_{2,1}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_{2,1}
%
%    f_2 is the indicator function of the set S define by Ax-b||_2 < epsilon
%     We define the prox of f_2 as 
%
%        prox_{f2,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma i_S( x ),
%
%     with i_S(x) is zero if x is in the set S and infinity otherwise.
%     This previous problem has an identical solution as:
%
%        argmin_{z} ||x - z||_2^2   s.t.  ||b - A z||_2 < epsilon
%
%     It is simply a projection on the B2-ball.
%     A is the measurement matrix (random Gaussian distribution)
%
%   The theoretical number of measurements M is computed with respect of the size of the
%   signal N and the sparsity level K:
%   
%      M=K*max(4,ceil(log(N)))   
%
%   Since we add some new information, we will try to reduce the number of
%   measurements by a factor p:
%
%      M=K*max(4/p,ceil(log(N)/p))   
%
%   With this number of measurements, we hope that the algorithm will
%   perform  a perfect reconstruction.
%
%   Results
%   -------
%
%   Figure 1: Results of the algorithm
%
%      This figure shows the original signal and the reconstruction done 
%      thanks to the algorithm and the measurements. The number of
%      measurements is M=900, the length of the signal N=5000, K=100, p=4. This is
%      equivalent to a compression ratio of 16.67. The elements are grouped
%      by 10.
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185--212, 2011.
%     
%     P. Combettes and J. Pesquet. A douglas--rachford splitting approach to
%     nonsmooth convex variational signal recovery. Selected Topics in Signal
%     Processing, IEEE Journal of, 1(4):564--574, 2007.
%     
%     F. Bach, R. Jenatton, J. Mairal, and G. Obozinski. Optimization with
%     sparsity-inducing penalties. arXiv preprint arXiv:1108.0775, 2011.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/demos/demo_compress_sensing3.php

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
% Date: november 2012


%% Initialisation

clear;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2; % verbosity level

%% Creation of the problem

g = 10;         % number of element per group
N = 5000;       % Size of the signal
K = g * 10;     % Sparsity level
p = 4;          % Gain with respect to the "traditional compression ratio"
R = max(ceil(4/p),ceil(log(N)/p));    % Constant 
fprintf('The compression ratio is: %g\n',N/(R*K));

% Mesurements matrix
A = randn(R * K, N);

% Create a K sparse signal
x = zeros(N, 1);
I2 = randperm(N/g)*g;
I = zeros(size(x));

for i=0:N/g-1
    I(i*g+1:(i+1)*(g)) = I2(i+1) * ones(g,1) - (1:g)' + 1;
end

x(I(1:K)) = randn(K, 1); % take a normal distribution (adapted to L12 norm)
x = x / norm(x);

% Create groups;
g_d = 1:N;
g_t = g*ones(1, N/g);

% Measurements
y = A * x;

%% Defining proximal operators

% Define the prox of f2 see the function proj_B2 for more help
operatorA = @(x) A * x;
operatorAt = @(x) A' * x;
epsilon2 = 1e-5;
param_proj.epsilon = epsilon2;
param_proj.A = operatorA;
param_proj.At = operatorAt;
param_proj.y = y;
param_proj.tight = 0;
param_proj.nu = norm(A)^2;
param_proj.verbose = verbose - 1;
f2.prox = @(x,T) proj_b2(x, T, param_proj);
f2.eval = @(x) norm(A*x - y)^2;

% setting the function f1
param_l21.verbose = verbose - 1;
param_l21.g_d = g_d;
param_l21.g_t = g_t;
f1.prox = @(x, T) prox_l21(x, T, param_l21);
f1.eval = @(x) norm_l21(x,g_d,g_t);   

%% solving the problem

% setting different parameters for the simulation
param_solver.verbose = verbose; % display parameter
param_solver.maxit = 300;       % maximum number of iterations
param_solver.tol = 1e-4;        % tolerance to stop iterating
param_solver.gamma = 1e-2;      % stepsize

% solving the problem
sol = solvep(zeros(N,1), {f1, f2}, param_solver);

%% displaying the result
figure;
plot(1:N, x, 'o', 1:N, sol, 'xr');
legend('Original signal', 'Reconstructed signal');
    
%% Closing the toolbox
close_unlocbox();

