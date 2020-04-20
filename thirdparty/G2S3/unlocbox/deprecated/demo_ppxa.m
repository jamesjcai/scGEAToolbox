%DEMO_PPXA  Example of use of the PPXA solver
%
%   We present an example of the ppxa solver through an image
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
%      This figure shows the original image(The cameraman). 
%
%   Figure 2: Depleted image
%
%      This figure shows the image after the application of the mask. Note
%      that 70% of the pixels have been removed.
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
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/deprecated/demo_ppxa.php

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

% Author: William Guicquero, Nathanael Perraudin
% Date: fev 23 2012
%

%% Initialisation

clear all;
close all;

% Loading toolbox
init_unlocbox();

verbose = 2;    % verbosity level

%% Defining the problem

% Original image
im_original = cameraman();

% Creating the problem
A = rand(size(im_original));
A = (A > 0.7);
% Depleted image
b=A.*im_original;
size_im=size(b);

  
%% Defining proximal operators

% setting the function f1 (norm TV)
param_tv.verbose = verbose - 1;
param_tv.maxit = 50;

f{1}.prox=@(x, T) reshape(prox_tv(reshape(x,size_im), T, param_tv),[],1);
f{1}.eval=@(x) norm_tv(reshape(x,size_im));   

% setting the function f2 (proj_b2)
operatorA=@(x) A.*x;
epsilon_ind=0e-2;
param_proj.epsilon=epsilon_ind;
param_proj.A=operatorA;
param_proj.At=operatorA;
param_proj.y=b;
param_proj.verbose = verbose - 1;

f{2}.prox=@(x,T) reshape(proj_b2(reshape(x,size_im),T,param_proj),[],1);
f{2}.eval=@(x) norm(A(:).*x(:)-b(:))^2;

%% Solving the problem

% setting different parameter for the simulation
param_solver.verbose = verbose;     % display parameter
param_solver.maxit = 200;           % maximum iteration
param_solver.tol = 10e-5;           % tolerance to stop iterating
param_solver.gamma = 0.5;            % stepsize

% solving the problem
sol = reshape(ppxa(b(:),f,param_solver),size(b));

%% displaying the result
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();


