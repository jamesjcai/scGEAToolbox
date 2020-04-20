% plot_some_images(X, pix_y, pix_x, n_rows, n_cols)
%
%
% example: 
% 
% load att_faces
% figure; plot_some_images(att_faces, 112, 92, 10, 10)
%
% load coil_20
% figure; plot_some_images(X, 128, 128, 15, 15)
%
% load coil_5_unprocessed
% figure; plot_some_images(X, 416, 448, 5, 5)
%
%
% code author: Vassilis Kalofolias
% date: 2015
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/plot_some_images.php

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

function plot_some_images(X, pix_y, pix_x, n_rows, n_cols)

if nargin<5
   n_cols = round(sqrt(size(X,2))*pix_x/pix_y);
end

if nargin<4
   n_rows = floor(size(X,2)/n_cols);
end

% the big image
I = zeros(pix_y * n_rows, pix_x * n_cols);

% index of image
k = 0;
for i = 1:n_rows
    for j = 1:n_cols
        k = k + 1;
        I(1+(i-1)*pix_y : i*pix_y, 1+(j-1)*pix_x : j*pix_x) = reshape(X(:, k), pix_y, pix_x);
    end
end
%imshow(I, []);
imagesc(I);
colormap gray


