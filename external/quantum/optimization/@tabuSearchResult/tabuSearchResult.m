classdef (Sealed) tabuSearchResult < quantum.internal.optimization.AbstractAlgorithmResult
%TABUSEARCHRESULT Tabu Search Result
% 
%   RES = TABUSEARCHRESULT creates an object to hold the result from
%   calling TABUSEARCH on a problem.

%   Copyright 2022-2023 The MathWorks, Inc.

    properties (SetAccess = private, GetAccess = public)
       AllX
       AllFunctionValues
       BestX
       BestFunctionValue
       TimeFirstBest
       TabuIterations
       NumTabuSearchCalls
       Message
       LongestStallTime
    end

    methods (Hidden)
        function obj = update(obj, xsol, fval, ~, output)
            if ~isempty(xsol)
                obj.AllX = output.AllX;
                obj.AllFunctionValues = output.AllFunctionValues;
                obj.BestX = xsol;
                obj.BestFunctionValue = fval;              
                obj.ClockTime = seconds(output.ElapsedTime);
                obj.TimeFirstBest = seconds(output.TimeFirstBest);
                obj.TabuIterations = output.TabuIterations;
                obj.NumTabuSearchCalls = output.NumTabuSearchCalls;
                obj.Message = output.Message;
                obj.LongestStallTime = seconds(output.LongestStallTime);
            end
        end
        function x = retrieveBestPoint(obj)
            x = obj.BestX;
        end  

        function obj = adjustForConstant(obj, constant)
            obj.BestFunctionValue = obj.BestFunctionValue + constant;
            obj.AllFunctionValues = obj.AllFunctionValues + constant; 
        end
    end
end