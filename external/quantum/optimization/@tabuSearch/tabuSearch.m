classdef (Sealed) tabuSearch < quantum.internal.optimization.tabusearch.MultiCallTabuSearch
%TABUSEARCH Tabu search
%
%   TABUSEARCH is an implementation of the Iterated Tabu Search algorithm
%   described in:
%
%   Iterated Tabu Search for the Unconstrained Binary Quadratic
%   Optimization Problem, Palubeckis, G., Informatica (2006), 17(2),
%   279-296

%   Copyright 2022-2023 The MathWorks, Inc.

    properties (Hidden)
        RestartCandidateSize (1, 1) {quantum.internal.optimization.tabusearch.MultiCallTabuSearch.mustBeDouble, mustBeInteger, mustBePositive} = 5 
        RestartLowerBound (1, 1) {quantum.internal.optimization.tabusearch.MultiCallTabuSearch.mustBeDouble, mustBeInteger, mustBePositive} = 10
        RestartVariableFraction (1, 1) {quantum.internal.optimization.tabusearch.MultiCallTabuSearch.mustBeDouble, mustBeNonnegative, mustBeReal, mustBeLessThanOrEqual(RestartVariableFraction, 1)} = 0.1
    end

    properties (Hidden, SetAccess = private, GetAccess = public)
        TabuSearchVersion = 1
    end

    methods

        function obj = tabuSearch(NameValueArgs)

            arguments 
                NameValueArgs.?tabuSearch
            end

            % Set any specified name-value pairs
            fNames = fieldnames(NameValueArgs);
            for i = 1:numel(fNames)
                obj.(fNames{i}) = NameValueArgs.(fNames{i});
            end

        end
	
        function result = doRun(obj, Q, constantForDisplay)
            
            [xsol, fval, exitflag, output] = search(obj, Q, constantForDisplay);
            result = tabuSearchResult;
            result = result.update(xsol, fval, exitflag, output);
            
        end
    end
    
end
