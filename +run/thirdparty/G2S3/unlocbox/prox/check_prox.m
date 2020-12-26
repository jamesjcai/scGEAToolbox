function is_ok = check_prox(f_eval, prox, size_input, x0, verbose)
% is_ok = check_gradients(f_eval, grad, size_input, x0, mode, verbose, symmetric)
%   f:              function to be checked
%   prox:           gradient of f
%   size_input:     if f_eval(X) is valid, then size_input = size(X);
%   x0:             the point at which the gradient is checked. If not
%                   provided, this is a randn matrix or vector. This is
%                   very useful when additional constraints are used for f,
%                   for example when it is only defined for positive
%                   values.
%   verbose:        0: no output, 1: details, 2: more details
%
%
% code author: Vassilis Kalofolias
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/prox/check_prox.php

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

%% TODO
if nargin < 4
    x0 = [];
end
if nargin < 5
    verbose = 1;
end

gamma_vals = logspace(-4, 3, 5);
n_prox = 10;

if isempty(x0)
    x0 = randn(size_input);
end

f = zeros(length(gamma_vals)+1, 1);
shrink_ok = zeros(size(gamma_vals));

for i_g = 1:length(gamma_vals)
    gamma = gamma_vals(i_g);
    fprintf('--Checking for gamma = %f\n', gamma);
    x_i = x0;
    f(1) = f_eval(x0);
    for i = 1:n_prox
        x_i = prox(x_i, gamma);
        f(i+1) = f_eval(x_i);
    end
    if verbose > 1
        disp(diff(f));
    end
    if all(diff(f) <= 0)
        shrink_ok(i_g) = 1;
    else
        shrink_ok(i_g) = 0;
    end
end

is_ok = all(shrink_ok == 1);


if verbose
    fprintf('Tested for gamma values: \n');
    disp(gamma_vals);
    disp(shrink_ok);
    %fprintf('Theoretical derivative: %f, approximation: %f\n', df, df_approx);
    if is_ok, fprintf('proximal operator seems OK\n'); end
end


