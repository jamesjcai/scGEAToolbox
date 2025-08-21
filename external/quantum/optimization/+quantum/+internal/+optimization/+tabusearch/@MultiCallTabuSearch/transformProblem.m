function [Cp, dp] = transformProblem(Cp, dp, x)
%TRANSFORMPROBLEM Transform problem
%
%   [CP, DP] = TRANSFORMPROBLEM(CP, DP, X) transforms the problem given by
%   the quadratic terms, CP and the linear terms, DP by mapping the current
%   point X to zero.
%
%   This is an implementation of equations (3) and (4) in:
%
%   Iterated Tabu Search for the Unconstrained Binary Quadratic
%   Optimization Problem, Palubeckis, G., Informatica (2006), 17(2),
%   279-296

%   Copyright 2022 The MathWorks, Inc.

% Transform linear term
colIdx = find(x == 1);
addTerm = 2*sumSelectedColumns(Cp, colIdx);
dp = (1 - 2*x).*(dp + addTerm);

% Transform quadratic term
Cp = calculateQuadraticTransform(Cp, x);

end
