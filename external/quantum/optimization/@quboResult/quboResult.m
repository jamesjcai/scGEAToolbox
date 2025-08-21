classdef (Sealed) quboResult
%QUBORESULT QUBO Result 
%
%   RES = QUBORESULT contains the best point and corresponding function value
%   based on the results available from individual algorithm. QUBOResult
%   provides all the information from the algorithm so no re-run is
%   necessary to obtain additional result data.
%   
%   Copyright 2022-2023 The MathWorks, Inc.     
    properties (SetAccess = private, GetAccess = public)
        BestX
        BestFunctionValue
        AlgorithmResult
    end
    methods
        function obj = quboResult(anAlgorithmResult, theQUBO)
            p = inputParser;
            validAlgorithmResult = @(x) isa(x,'quantum.internal.optimization.AbstractAlgorithmResult');
            addRequired(p,'anAlgorithmResult',validAlgorithmResult);
            validQUBO = @(x) isa(x,'qubo');
            addRequired(p,'theQUBO',validQUBO);            
            parse(p,anAlgorithmResult, theQUBO);

            obj.AlgorithmResult = p.Results.anAlgorithmResult;
            obj.BestX = obj.AlgorithmResult.retrieveBestPoint;
            if isempty(obj.BestX)
                obj.BestFunctionValue = [];
            else
                obj.BestFunctionValue = p.Results.theQUBO.evaluateObjective(obj.BestX);
            end
        end
    end
end