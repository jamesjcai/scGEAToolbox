function [bestx, bestf, Cp, dp, x, iter] = callTabuOnTransformedProblem(obj, Cp, dp, x, initf, bestx, bestf, startTime)
%CALLTABUONTRANSFORMEDPROBLEM Search for an improved solution using Tabu Search
%
%   [BESTX, BESTF, CP, DP, X, ITER] = CALLTABUONTRANSFORMEDPROBLEM(OBJ, CP, DP, X, INITF,
%   BESTX, BESTF, STARTTIME) calls Tabu Search on the transformed quadratic
%   unconstrained binary optimization (QUBO) problem, CP, DP.
%
%   This is an implementation of algorithm TS in the following reference
%
%   Iterated Tabu Search for the Unconstrained Binary Quadratic
%   Optimization Problem, Palubeckis, G., Informatica (2006), 17(2),
%   279-296

%   Copyright 2022-2023 The MathWorks, Inc.

% Number of variables
numVar = numel(x);

% Initialise tabu values
L = zeros(numVar, 1);

% Set tabu tenure
TabuTenure = ceil(min(20, numVar/4));

% Initialize run time defaults
[obj, maxStallIter] = initializeRunTimeDefaults(obj, numVar);

% Initialize iteration counters
iter = 0;
numStallIter = 0;

% Initialize fbest
fbar = initf;

% Main loop
done = false;
while ~done

    % Examine all non-tabu variables
    % This is the neigbourhood search phase
    [idxBestVar, doLocalSearch, nbhdIter] = quantum.internal.optimization.tabusearch.MultiCallTabuSearch.neigborhoodsearch(L, dp, numVar, fbar, bestf);
    numStallIter = numStallIter + nbhdIter;
    iter = iter + nbhdIter;

    % Flip idxBestVar-th variable 
    x(idxBestVar) = 1 - x(idxBestVar);
    fbar = fbar + dp(idxBestVar);
    [Cp, dp] = quantum.internal.optimization.tabusearch.MultiCallTabuSearch.updatemodel(Cp, dp, idxBestVar);

    % Apply the local search procedure
    if doLocalSearch
        numStallIter = 0;
        [localx, flocal, localiter, Cp, dp] = quantum.internal.optimization.tabusearch.MultiCallTabuSearch.localsearch(x, Cp, dp);
        iter = iter + localiter;
        if flocal > 0
            % Record variable associated with the move
            x = localx;
            fbar = fbar + flocal;
        end
        bestf = fbar;
        bestx = x;
    end

    % Reduce all tabu values by one
    L = L - 1;
    L(L < 0) = 0;

    % Tabu the chosen variable
    L(idxBestVar) = TabuTenure;

    % Check whether it's time to stop
    done = iter > obj.MaxIterations || toc(startTime) > obj.MaxTime || ...
        numStallIter > maxStallIter;

end

end






