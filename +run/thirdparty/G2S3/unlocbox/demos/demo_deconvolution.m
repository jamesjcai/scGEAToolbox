%DEMO_DECONVOLUTION Deconvolution demonstration (Debluring)
%
%   Here we try to deblur an image through a deconvolution problem. The
%   convolution operator is the blur
%   The problem can be expressed as this
%
%        argmin  ||Ax-b||^2 + tau*||H(x)||_1
%
%   Where b is the degraded image, I the identity and A an operator representing the blur.
%
%   H is a linear operator projecting the signal in a sparse
%   representation. Here we worked with wavelet. 
%
%   Warning! Note that this demo require the LTFAT toolbox to work.
%
%   We set 
%
%    f_1(x)=||H(x)||_{1}
%     We define the prox of f_1 as: 
%
%        prox_{f1,gamma} (z) = argmin_{x} 1/2 ||x-z||_2^2  +  gamma  ||H(z)||_1
%
%    f_2(x)=||Ax-b||_2^2
%     We define the gradient as: 
%
%        grad_f(x) = 2 A^*(Ax-b)
%
%   Results
%   -------
%
%   Figure 1: Original image
%
%      This figure shows the original lena image. 
%
%   Figure 2: Depleted image
%
%      This figure shows the image after the application of the blur.
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
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/demos/demo_deconvolution.php

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
% Date: sept 30 2013
%

%% Initialisation

clear;
close all;

init_unlocbox;

verbose = 2;    % verbosity level
clim = [0 1];   % limits for colors

%% Define the problem

% Original image
im_original = cameraman();

% Create a blur
sigma = 0.1;
[x, y] = meshgrid(linspace(-1, 1, length(im_original)));
r = x.^2 + y.^2;
G = exp(-r/(2*sigma^2));

A=@(x) real(ifft2(fftshift(G).*(fft2(x))));
At=@(x) real(ifft2(conj(fftshift(G)).*(fft2(x))));

% Depleted image
b=A(im_original);

%% Define proximity operators

% to deblur with wavelet
tau = 0.0003; %parameter for the problem


% Use LTFAT
L=8;
W = @(x) fwt2(x,'db1',L);
Wt = @(x) ifwt2(x,'db1',L);

param_l1.verbose = verbose-1;
param_l1.tight = 1;
param_l1.At = Wt;
param_l1.A = W;

f.prox=@(x, T) prox_l1(x, T*tau, param_l1);
f.eval=@(x) tau*sum(sum(abs(W(x))));   

% % to deblur with TV
% tau = 0.001;
% param_tv.verbose = verbose-1;
% f.prox=@(x, T) prox_tv(x, T*tau, param_tv);
% f.eval=@(x) tau*tv_norm(x);   

param_proj.maxit = 10;
param_proj.epsilon = 0;
param_proj.tight = 0;
param_proj.nu = 2;
param_proj.A = A;
param_proj.At = At;
param_proj.y = b;
param_proj.verbose = verbose - 1;
f2.eval = @(x) norm(A(x) - b).^2;
f2.prox = @(x,T) proj_b2(x, T, param_proj);


%% Solve the problem

% setting different parameter for the simulation
param_solver.verbose=verbose; % display parameter
param_solver.maxit=300; % maximum iteration
param_solver.tol=10e-9; % tolerance to stop iterating
fig = figure(100);
param_solver.do_sol=@(x) plot_image(x,fig); % plotting plugin
% solving the problem
sol=solvep(b, {f, f2}, param_solver);
close(fig)
%% displaying the result
imagesc_gray(im_original, 1, 'Original image',111,clim) 
imagesc_gray(b, 2, 'Depleted image',111,clim) 
imagesc_gray(sol, 3, 'Reconstructed image',111,clim) 
    

%% Closing the toolbox
close_unlocbox();



