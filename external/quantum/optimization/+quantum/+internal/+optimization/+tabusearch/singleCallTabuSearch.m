function [bestx, bestf] = singleCallTabuSearch(Q, numIterations, x0)
%SINGLECALLTABUSEARCH Search for a solution to a QUBO with a single call to
%Tabu search

%   [BESTX, BESTF] = SINGLECALLTABUSEARCH(C, NUMITERATION, X0) calls Tabu
%   searh on the Quadratic Unconstrained Binary problem defined by quadratic
%   term Q and linear term d. 

%   This is an implementation of algorithm TS in the following reference
%
%   Iterated Tabu Search for the Unconstrained Binary Quadratic
%   Optimization Problem, Palubeckis, G., Informatica (2006), 17(2),
%   279-296

%   Copyright 2022-2023 The MathWorks, Inc.

% Create TabuSearch object
obj = tabuSearch;
obj.MaxIterations = numIterations; 

% QuboProblem minimizes, whereas TabuSearch maximizes
Q = -Q;

% Get initial point
if nargin == 3
    x = x0(:);
else
    x = obj.initialPoint(size(Q,1));
end

% Prepare and transform the QUBO
[Cp,dp,x,bestx,bestf] = obj.prepareQuboForTabu(Q, x);

[bestx, bestf] = obj.callTabuOnTransformedProblem(Cp, dp, x, bestf, bestx, bestf, tic);

bestf = -bestf;

end