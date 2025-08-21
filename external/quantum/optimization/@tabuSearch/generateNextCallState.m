function [Cp, dp, x] = generateNextCallState(obj, Cp, dp, bestx)
%GENERATENEXTCALLSTATE Generate state for next call to Tabu Search
%
%   [CP, DP, X] = GENERATENEXTCALLSTATE(OBJ, CP, DP, CURRX) is a class to
%   generate state for the next call to TabuSearch using the Get Start
%   Point (GSP) algorithm in
%
%   Iterated Tabu Search for the Unconstrained Binary Quadratic
%   Optimization Problem, Palubeckis, G., Informatica (2006), 17(2),
%   279-296

%   Copyright 2022 The MathWorks, Inc.

numVar = numel(bestx);
restartSize = min(obj.RestartCandidateSize, numVar);
numRestartVars = ceil(numVar*obj.RestartVariableFraction);
rbarLower = min(obj.RestartLowerBound, numRestartVars);
rbarUpper = numRestartVars;
rbar = randi([rbarLower, rbarUpper]);
[x, Cp, dp] = getstartpoint(bestx, Cp, dp, restartSize, rbar);

end

function [x, Cp, dp] = getstartpoint(x, Cp, dp, b, rbar)

% Step 1: Initialize
r = 0;
numVar = numel(x);
I = 1:numVar;
idxRemove = false(numVar, 1);

while r < rbar

    % Instead of removing thisIdx from I, set dp(thisIdx) to -inf in this
    % routine. This has the same effect as removing thisIdx from I. Also
    % means that we don't have to map indices.
    dpSort = dp;
    dpSort(idxRemove) = -inf;

    % Step 2: Pick the b largest coefficients in dp
    [~, idx] = sort(dpSort, 'descend');

    % Select variable to flip
    q = randi([1 b]);
    thisIdx = I(idx(q));

    % Step 3:  Flip selected variable. Update model.
    x(thisIdx) = 1 - x(thisIdx);
    [Cp, dp] = quantum.internal.optimization.tabusearch.MultiCallTabuSearch.updatemodel(Cp, dp, thisIdx);
    idxRemove(thisIdx) = true;

    % Step 4: Increment r by 1
    r = r + 1;

end

end