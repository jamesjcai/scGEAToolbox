function [x, flocal, iter, Cp, dp] = localsearch(x, Cp, dp)
%LOCALSEARCH Find a local optimum 
%
%   [X, FLOCAL, NUMITER, CP, DP] = LOCALSEARCH(X, CP, DP) finds a local
%   optimum from the start point X, using a simple ascent procedure. The
%   current state of the model is passed in and out via CP and DP. The
%   local optimum is returned in X, the improvement in the function value
%   in FLOCAL and the number of iterations taken in NUMITER.
% 
%   NOTE: This is an implementation of the LS procedure in:
%
%   Iterated Tabu Search for the Unconstrained Binary Quadratic
%   Optimization Problem, Palubeckis, G., Informatica (2006), 17(2),
%   279-296

%   Copyright 2022 The MathWorks, Inc.

flocal = 0;
done = false;
iter = 0;
while ~done

    done = true;
    for k = 1:numel(x)
        iter = iter + 1;
        if dp(k) > 0
            done = false;
            x(k) = 1 - x(k);
            flocal = flocal + dp(k);
            [Cp, dp] = quantum.internal.optimization.tabusearch.MultiCallTabuSearch.updatemodel(Cp, dp, k);
        end
    end


end
