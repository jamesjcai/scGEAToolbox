%DEMO_DOUGLAS_RACHFORD  Example of use of the douglas_rachford solver
%
%   We present an example of the douglas_rachford solver through an image
%   reconstruction problem.
%   The problem can be expressed as this
%
%       argmin ||x||_TV s.t ||b-Ax||_2 < epsilon
%
%   Where b is the degraded image, I the identity and A an operator representing the mask.
%
%   Note that the constraint can be inserted in the objective function
%   thanks to the help of the indicative function. Then we recover the
%   general formulation used for the solver of this toolbox.
%
%   We set 
%
%    f_1(x)=||x||_{TV}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_TV
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
%
%   Results
%   -------
%
%   Figure 1: Original image
%
%      This figure shows the original Lena image. 
%
%   Figure 2: Depleted image
%
%      This figure shows the image after the application of the mask. Note
%      that 85% of the pixels have been removed.
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
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/demos/demo_douglas_rachford.php

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

clear;
close all;

init_unlocbox;

verbose = 2;    % Verbosity level

%% Creation of the problem

% Original image
im_original = lena();   

% Creating the problem
A = rand(size(im_original));
A = A > 0.85;

% Depleted image
b = A .* im_original;

%% Defining proximal operators

% Define the prox of f2 see the function proj_B2 for more help
operatorA = @(x) A.*x;
epsilon2 = 0;
param_proj.epsilon = epsilon2;
param_proj.A = operatorA;
param_proj.At = operatorA;
param_proj.y = b;
param_proj.verbose = verbose - 1;
f2.prox=@(x,T) proj_b2(x, T, param_proj);
f2.eval=@(x) eps;

f2.prox = @(x,T) (x - A.*x) + A.*b;


% setting the function f1 (norm TV)
param_tv.verbose = verbose - 1;
param_tv.maxit = 50;

f1.prox = @(x, T) prox_tv(x, T, param_tv);
f1.eval = @(x) norm_tv(x);   

%% solving the problem

% setting different parameter for the simulation
paramsolver.verbose = verbose;  % display parameter
paramsolver.maxit = 100;        % maximum iteration
paramsolver.tol = 10e-7;        % tolerance to stop iterating
paramsolver.gamma = 0.1;        % stepsize

% To see the evolution of the reconstruction
fig = figure(100);
paramsolver.do_sol = @(x) plot_image(x,fig);

sol = douglas_rachford(b,f1,f2,paramsolver);

close(100);

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();


