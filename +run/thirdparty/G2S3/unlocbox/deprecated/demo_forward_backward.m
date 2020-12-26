%DEMO_FORWARD_BACKWARD  Example of use of the forward_backward solver
%
%   We present an example of the forward_backward solver through an image
%   reconstruction problem.
%   The problem can be expressed as this
%
%        argmin ||Ax-b||^2 + tau*||x||_TV
%
%   Where b is the degraded image, I the identity and A an operator representing the mask.
%
%   We set 
%
%    f_1(x)=||x||_{TV}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_TV
%
%    f_2(x)=||Ax-b||_2^2
%     We define the gradient as: 
%
%        grad_f(x) = 2 * A^*(Ax-b)
%
%
%   Results
%   -------
%
%   Figure 1: Original image
%
%      This figure shows the original image (The cameraman).  
%
%   Figure 2: Depleted image
%
%      This figure shows the image after the application of the mask. Note
%      that 50% of the pixels have been removed.
%
%   Figure 3: Reconstructed image
%
%      This figure shows the reconstructed image thanks to the algorithm.
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185--212, 2011.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/deprecated/demo_forward_backward.php

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
% Date: November 2012
%

%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

tau = 0.1; % Regularization parameter

% Original image
im_original = cameraman();

% Creating the problem
A = rand(size(im_original));
A = (A > 0.5);

% Depleted image
b = A .* im_original;

%% Defining proximal operators

% setting the function f2 
f2.grad=@(x) 2 * A .* (A.*x - b);
f2.eval=@(x) norm(A(:).*x(:)-b(:))^2;

% setting the function f1 (norm TV)
param_tv.verbose=verbose - 1;
param_tv.maxit=50;

f1.prox=@(x, T) prox_tv(x, T*tau, param_tv);
f1.eval=@(x) tau * norm_tv(x);   

%% Solving the problem

% setting different parameter for the simulation
param_solver.verbose = verbose; % display parameter
param_solver.maxit = 50;        % maximum iteration
param_solver.gamma = 0.5;       % stepsize (beta is equal to 2)
param_solver.tol = 1e-6;        % Tolerance to stop iterating

% solving the problem
t1 = tic;
sol = forward_backward(b,f1,f2,param_solver);
time1 = toc(t1);


%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');  

%% Closing the toolbox
close_unlocbox();

