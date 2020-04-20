function fg = add_smooth_functions(Fg)
%ADD_SMOOTH_FUNCTIONS sum of smooth function
%   Usage: fg = add_smooth_functions(Fg);
%    
%   Input parameters:
%       Fg  : cell array of function (cell array of struct)
%   Output parameters:
%       fg  : function (struct)
%
%   This function takes a cell array of smooth functions and transform into
%   a single smooth function. The array can be empty.
%
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/utils/add_smooth_functions.php

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


if numel(Fg)==0
    fg.grad = @(x) 0;
    fg.eval = @(x) 0;
    fg.beta = 0;
elseif numel(Fg) == 1
    fg = Fg{1};
else
    beta = 0;
    for ii = 1:numel(Fg)
        beta = beta + Fg{ii}.beta;
    end

    fg.grad = @(x) sum_grad(Fg,x);
    fg.eval = @(x) sum_eval(Fg,x);
    fg.beta = beta;
end

if ~(numel(fg.beta)==1)
    error('f.beta must be a scalar!')
end

end



function curr_grad = sum_grad(Fg,x)

    curr_grad = 0;
    for ii = 1:numel(Fg)
        curr_grad = curr_grad + Fg{ii}.grad(x);
    end
    
end

function curr_eval = sum_eval(Fg,x)

    curr_eval = 0;
    for ii = 1:numel(Fg)
        curr_eval = curr_eval + Fg{ii}.eval(x);
    end

end
