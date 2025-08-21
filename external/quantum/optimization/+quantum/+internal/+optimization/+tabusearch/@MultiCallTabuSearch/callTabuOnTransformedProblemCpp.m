function [bestx, bestf, Cp, dp, currx, iter] = callTabuOnTransformedProblemCpp(obj, Cp, dp, x, initf, bestx, bestf, fixedStruct, elapsedTime)
%CALLTABUONTRANSFORMEDPROBLEMCPP Search for an improved solution using Tabu Search
%
%   [BESTX, BESTF, CP, DP, X, ITER] = CALLTABUONTRANSFORMEDPROBLEM(OBJ, CP,
%   DP, X, INITF, BESTX, BESTF, FIXEDSTRUCT, ELAPSEDTIME) calls Tabu Search
%   on the transformed quadratic unconstrained binary optimization (QUBO)
%   problem, CP, DP.
%
%   This is an implementation of algorithm TS in the following reference
%
%   Iterated Tabu Search for the Unconstrained Binary Quadratic
%   Optimization Problem, Palubeckis, G., Informatica (2006), 17(2),
%   279-296

%   Copyright 2023 The MathWorks, Inc.

% Use workaround to add support package Windows DLLs to system path
quantum.internal.aws.addDLLLocationToSystemPathIfNecessary();

% Best point
bestpt.X = bestx;
bestpt.Fval = bestf;

% Inputs that vary between each call
varStruct.Values = 2*Cp.Values(:);
varStruct.LinearTerm = dp;
varStruct.X = x;
varStruct.Fval = initf;

% Initialize run time defaults
[obj, maxStallIter] = initializeRunTimeDefaults(obj, fixedStruct.NumVariables);

% Simple tabu search options
options.MaxIterations = obj.MaxIterations;
options.MaxStallIterations = maxStallIter;
options.MaxTime = obj.MaxTime - elapsedTime;

% Call the driver
[results, varStruct] = quantum.internal.optimization.tabusearch.mx_simpletabusearch(fixedStruct, varStruct, bestpt, options);

% Map outputs back 
bestx = results.X; 
bestf = results.Fval;
iter = results.Iterations;

% Map variable inputs back
Cp.Values = 0.5*varStruct.Values;
if strcmp(fixedStruct.Type, "dense")
    Cp.Values = reshape(Cp.Values, fixedStruct.NumVariables, fixedStruct.NumVariables);
end
dp = varStruct.LinearTerm;
currx = varStruct.X;






