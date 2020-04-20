%DEMO_PROX_MULTI_FUNCTIONS Demonstration of the proximal operator of a sum of function
%   
%   In this example we solve a image denoising problem. We remember the
%   reader that solving a proximal operator can be view like denoising a
%   signal. The function used is prox_sumg which compute the proximal
%   operator of a sum of function.
%
%   The problem can be expressed as this
%
%        argmin ||x-y||^2 + tau1*||x||_TV + tau2 * ||H(x)||_1
%
%   Where z is the degraded image.
%
%   H is a linear operator projecting the signal in a sparse
%   representation. Here we worked with wavelet. 
%
%   Warning! Note that this demo require the rwt(RICE WAVELET TOOLBOX) to work.
%
%   We set 
%
%    f_1(x)=||x||_{TV}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma ||z||_TV
%
%    f_2(x)=||H(x)||_{1}
%     We define the prox of f_1 as: 
%
%        prox_{f2,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||H(z)||_1
%
%
%   Results
%   -------
%
%   Figure 1: Original image
%
%      This figure shows the original cameraman image. 
%
%   Figure 2: Depleted image
%
%      This figure shows the image after the addition of noise.
%
%   Figure 3: Reconstruted image
%
%      This figure shows the denoised image thanks to the algorithm.
%
%   References:
%     P. Combettes and J. Pesquet. Proximal splitting methods in signal
%     processing. Fixed-Point Algorithms for Inverse Problems in Science and
%     Engineering, pages 185--212, 2011.
%     
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/deprecated/demo_prox_multi_functions.php

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

 
% Author: Nathanael Perraudin, Gilles Puy
% Date: sept 30 2011
%


%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

im_original = cameraman();

% Depleted image
sigma = 0.75;
b = im_original + sigma^2 * rand(size(im_original));

%% Defining proximal operators

% setting different parameter for the simulation
tau1 = 0.15; % regularization parameter for the TV norm
tau2 = 0.05; % regularization parameter for the wavelet

% setting the function f1 (TV Norm)
param_tv.verbose = verbose - 1;
param_tv.maxit = 100;
g1.prox=@(x, T) prox_tv(x, T*tau1, param_tv);
g1.eval=@(x) tau1 * norm_tv(x);   

% setting the function f2 (Wavelet)
L=8;
% h = daubcqf(2);
% H = @(x) mdwt(x,h,L);
% Ht = @(x) midwt(x,h,L);
H = @(x) fwt2(x, 'db1', L);
Ht = @(x) ifwt2(x, 'db1', L);
param_l1.verbose = verbose - 1;
param_l1.tight = 1;
param_l1.At = Ht;
param_l1.A = H;
g2.prox=@(x, T) prox_l1(x, T*tau2, param_l1);
g2.eval=@(x) tau2 * norm(x(:),1);  

%% Solving the problem

% Parameter for the prox
G={g2,g1};
param_sumg.G=G;
param_sumg.maxit = 100;
param_sumg.verbose = verbose;
param_sumg.tol = 1e-5;


% solving the problem
sol=prox_sumg(b,1,param_sumg);

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();

