function curr_norm = eval_function(fg,Fp,x,s,param)
%EVAL_FUNCTION internal evaluation function
%
%   Url: https://epfl-lts2.github.io/unlocbox-html/doc/solver/misc/eval_function.php

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

% 
% if nargin<5
%     param = struct;
%     param.fast_eval = 0;
% end

curr_norm = fg.eval(x); 


for ii = 1:length(Fp)
    
%     % Here we try to accelerate the evaluation of the function by using the
%     % results of the proximal operator. The proximal functions of the
%     % UNLocBoX returns two arguments. In the second one, final_eval can be
%     % used as an approximation of the norm.
%     if param.fast_eval ...
%       && isfield(s,'x_n') ...
%       && iscell(s.x_n) ...
%       && length(s.x_n) >= ii ...
%       && iscell(s.x_n{ii}) ...
%       && length(s.x_n{ii})>1 ...
%       && isfield(s.x_n{ii}{2},'final_eval')
%         curr_norm = curr_norm + s.x_n{ii}{2}.final_eval/param.gamma;
%     else 
        curr_norm = curr_norm + Fp{ii}.eval(x); 
%     end
end

end
