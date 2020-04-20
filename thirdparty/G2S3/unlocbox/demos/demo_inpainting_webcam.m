
% Author : Nathanaël Perraudin
%
% This demo capture a image with the webcam and perform some inpainting
% with it.
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/demos/demo_inpainting_webcam.php

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

clear;

close all;

init_unlocbox;

%% parameter

scale = 0.5; % scale for the image
imgcolor = 1;
p = 0.8 ;

verbose = 2;    % Verbosity level
maxit = 30;

%%

% Close previous cam if3 any
clear('cam')

camList = webcamlist;
% Is there any cam?

if numel(camList)>1
    fprintf('Please choose a webcam:\n');
    for ii = 1:numel(camList)
        fprintf(['  ',num2str(ii),') ',camList{ii},'\n'])
    end
    fprintf('Camera number: ')
    prompt = 1;
    numCam = str2num(input('','s'));
else
    numCam = 1;
end

%% Start the webcam

cam = webcam(numCam);
%% Aquire image

preview(cam);

fprintf('Push a button to aquire image... ')
pause;

% Acquire a single image.
rgbImage = snapshot(cam);
fprintf('   Done! \n')

rgbImage = imresize(rgbImage, scale);
if imgcolor
    im_original = double(rgbImage)/256;
else
    % Convert RGB to grayscale and rescaling.
    grayImage = double(rgb2gray(rgbImage))/256;
    % Original image
    im_original = grayImage; 
end

closePreview(cam)
clear('cam')

%%


%% Creation of the problem

 

% Creating the problem
A = rand(size(im_original,1),size(im_original,2));
A = A > p;

if imgcolor
   A = cat(3,A,A,A); 
end

% Depleted image
b = A .* im_original;


imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');

pause;

%% Defining proximal operators

% Define the prox of f2 see the function proj_B2 for more help
operatorA = @(x) A.*x;
operatorAt = @(x) A.*x;
epsilon2 = 0;
param_proj.epsilon = epsilon2;
param_proj.A = operatorA;
param_proj.At = operatorAt;
param_proj.y = b;
param_proj.verbose = verbose - 1;
f2.prox=@(x,T) proj_b2(x, T, param_proj);
f2.eval=@(x) eps;

f2.prox = @(x,T) (x - A.*x) + A.*b;


% setting the function f1 (norm TV)
param_tv.verbose = verbose - 1;
param_tv.maxit = maxit;

f1.prox = @(x, T) prox_tv(x, T, param_tv);
f1.eval = @(x) sum(norm_tv(x));   

%% solving the problem

% setting different parameter for the simulation
paramsolver.verbose = verbose;  % display parameter
paramsolver.maxit = maxit;        % maximum iteration
paramsolver.tol = 10e-7;        % tolerance to stop iterating
paramsolver.gamma = 0.1;        % stepsize

% To see the evolution of the reconstruction
fig = figure(100);
paramsolver.do_sol = @(x) plot_image(x,fig,1);

sol = douglas_rachford(b,f1,f2,paramsolver);

close(100);

%% displaying the result
sol(sol<0) = 0;
sol(sol>1) = 1;
imagesc_gray(im_original, 1, 'Original image');
imagesc_gray(b, 2, 'Depleted image');
imagesc_gray(sol, 3, 'Reconstructed image');
    

%% Closing the toolbox
close_unlocbox();

