classdef MultiCallTabuSearch < quantum.internal.optimization.AbstractAlgorithm
%MULTICALLTABUSEARCH Abstract class for multi-call Tabu Search algorithms 
%
%   MULTICALLTABUSEARCH is an abstract class to search for 
%   solution to a QUBO by calling a TabuSearch algorithm several times.

%   Copyright 2022-2023 The MathWorks, Inc.

    properties
        Display (1, 1) string {mustBeMember(Display, ["off", "final", "iter"])} = "off"
        MaxIterations {mustBeScalarOrEmpty, quantum.internal.optimization.tabusearch.MultiCallTabuSearch.mustBeDouble, mustBeInteger, mustBePositive} 
        MaxStallIterations {mustBeScalarOrEmpty, quantum.internal.optimization.tabusearch.MultiCallTabuSearch.mustBeDouble, mustBeInteger, mustBePositive} 
        MaxStallTime {mustBeScalarOrEmpty, quantum.internal.optimization.tabusearch.MultiCallTabuSearch.mustBeDouble, mustBePositive, mustBeReal, mustBeFinite} 
        MaxTime (1, 1) {quantum.internal.optimization.tabusearch.MultiCallTabuSearch.mustBeDouble, mustBePositive, mustBeReal, mustBeFinite} = 20
    end

    properties (Hidden)
        ShowIterativeDisplay = false
        ShowDiagnostics = false
        StoreQSparse = true
        Mu = 10000
        UseCppTabuEngine = true
        MaxSingleTabuCalls = Inf
        ReturnTransformedProblem = false
    end

    properties (Hidden, SetAccess = private, GetAccess = public)
        MultiCallTabuSearchVersion = 1
    end
    
    methods (Abstract)
        result = doRun(obj, quboProblem)
        [Cp, dp, x] = generateNextCallState(obj, Cp, dp, currx, varargin)
    end

    methods
        [Cp,dp,x,bestx,bestf] = prepareQuboForTabu(obj, C, x0)
        [bestx, bestf, Cp, dp, x, iter] = callTabuOnTransformedProblem(obj, Cp, dp, x, initf, bestx, bestf, startTime)
        [bestx, bestf, Cp, dp, x, iter] = callTabuOnTransformedProblemCpp(obj, Cp, dp, x, initf, bestx, bestf, startTime, fixedStruct)
        [bestx, bestf, exitflag, output] = search(obj, C, constantForDisplay)
        updateIterativeDisplay(obj, numTabuCalls, bestFval, elapsedTime, tabuIter)
        [obj, maxStallIter] = initializeRunTimeDefaults(obj, numVar)
    end

    methods (Static)
        [C, d] = moveDiagonal2Linear(C, d)
        [Cp, dp] = transformProblem(Cp, dp, x)
        [Cp, dp] = updatemodel(Cp, dp, idx)
        [idxBestVar, doLocalSearch, iter] = neigborhoodsearch(L, linTerm, numVar, fbar, bestf)
        [x, flocal, numIter, Cp, dp] = localsearch(x, Cp, dp)
        x0 = initialPoint(numVar);
        mustBeDouble(input);
    end

end
