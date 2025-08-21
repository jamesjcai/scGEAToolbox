classdef (Sealed) qaoaResult < quantum.internal.optimization.AbstractAlgorithmResult
%QAOARESULT Quantum Approximate Optimization Algorithm Result
% 
%   RES = QAOARESULT creates an object to hold the result from calling QAOA
%   on a problem.

%   Copyright 2024 The MathWorks, Inc.

    properties (SetAccess = private, GetAccess = public)
       BestX
       BestFunctionValue
       AllX
       AllFunctionValues
       AlgorithmInformation = struct("AllProbabilities", [], "BestCircuitAngles", [], "BestExpectedValue", [])
    end

    methods (Hidden)
        function obj = update(obj, xsol, fval, exitflag, output)
            if ~isempty(xsol)
                obj.BestX = xsol;
                obj.BestFunctionValue = fval;
                obj.AllX = output.AllX;
                obj.AllFunctionValues = output.AllFunctionValues;
                obj.AlgorithmInformation.AllProbabilities = output.AllProbabilities;
                obj.AlgorithmInformation.BestCircuitAngles = output.BestCircuitAngles;
                obj.AlgorithmInformation.BestExpectedValue = output.BestExpectedValue;
                obj.AlgorithmInformation.Exitflag = exitflag;
                obj.AlgorithmInformation.Output = output.ClassicalOutput;
                obj.AlgorithmInformation.Circuit = output.Circuit;
                obj.ClockTime = output.ClockTime;
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