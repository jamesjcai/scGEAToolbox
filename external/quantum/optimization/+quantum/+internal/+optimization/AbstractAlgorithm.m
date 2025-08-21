classdef AbstractAlgorithm < matlab.mixin.Heterogeneous
%ABSTRACTALGORITHM Base class for Algorithm classes

%   Copyright 2022-2023 The MathWorks, Inc.

    methods(Abstract=true)
        result = doRun()
    end

    methods
        function result = run(obj, quboProblem)
            [Q, constant] = unwrapQubo(obj, quboProblem);
            result = doRun(obj, Q, constant);
            result = result.adjustForConstant(constant);
        end

        function [Q, constant] = unwrapQubo(~, quboProblem)

            % Add the linear term to the diagonal
            if issparse(quboProblem.QuadraticTerm)
                N = quboProblem.NumVariables;
                Q = quboProblem.QuadraticTerm + spdiags(quboProblem.LinearTerm, 0, N, N);
            else
                Q = quboProblem.QuadraticTerm + diag(quboProblem.LinearTerm);
            end

            % Extract constant
            constant = quboProblem.ConstantTerm;

        end
    end
end