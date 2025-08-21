classdef (Sealed) qaoa < quantum.internal.optimization.AbstractAlgorithm & matlab.mixin.CustomDisplay
%QAOA Quantum Approximate Optimization Algorithm
%
%   QAOA is an implementation of the Quantum Approximate Optimization
%   Algorithm described in:
%
%   A Quantum Approximate Optimization Algorithm, Edward Farhi and Jeffrey
%   Goldstone and Sam Gutmann, 2014, https://arxiv.org/abs/1411.4028 

%   Copyright 2024 The MathWorks, Inc.

    properties
        InitialAngles {mustBeNumeric, mustBeFinite, mustBeReal} = []
        NumLayers (1,1) {mustBeNumeric, mustBeScalarOrEmpty, mustBeFinite, mustBePositive, mustBeInteger} = 1
        NumShots (1,1) {mustBeNumeric, mustBeScalarOrEmpty, mustBeFinite, mustBePositive, mustBeInteger} = 100         
    end

    properties (Dependent)
        OptimizationSolver {mustBeMember(OptimizationSolver, ["fminsearch", "fmincon", "surrogateopt"])} 
        OptimizationSolverOptions {mustBeValidOptimoptionsOrStructOrEmpty}
    end

    properties (Hidden, SetAccess = private)
        OptimizationConfiguration = struct("Solver", "fminsearch", "Options", [])
    end

    properties (Hidden, Constant)
        ValidSolvers = ["surrogateopt", "fmincon", "fminsearch"]
    end

    methods

        function obj = qaoa(NameValueArgs)

            arguments 
                NameValueArgs.?qaoa
            end

            % Set any specified name-value pairs
            fNames = fieldnames(NameValueArgs);
            for i = 1:numel(fNames)
                obj.(fNames{i}) = NameValueArgs.(fNames{i});
            end

        end

        function result = doRun(obj, Q, constantForDisplay)
            
            [xsol, fval, exitflag, output] = search(obj, Q, constantForDisplay);
            result = qaoaResult;
            result = result.update(xsol, fval, exitflag, output);
            
        end

    end

    % Set/get methods for dependent properties
    methods 
        function solver = get.OptimizationSolver(obj)            
            solver = obj.OptimizationConfiguration.Solver;
        end

        function obj = set.OptimizationSolver(obj, solver)
            obj.OptimizationConfiguration.Solver = solver;
            obj.OptimizationConfiguration.Options = [];
        end

        function options = get.OptimizationSolverOptions(obj)
            options = obj.OptimizationConfiguration.Options;
        end

        function obj = set.OptimizationSolverOptions(obj, options)
            
            % Check to see if the supplied OptimizationSolverOptions are
            % valid for the current OptimizationSolver. 
            % If the solver is fminsearch, the options must be a struct.
            % Otherwise the options must be an optimoptions object for the
            % given solver.
            isOptionsForSolver = ...
                ( isstruct(options) && strcmp(obj.OptimizationConfiguration.Solver, "fminsearch") ) || ...
                strcmp(obj.OptimizationConfiguration.Solver, options.SolverName);
            if isOptionsForSolver
                obj.OptimizationConfiguration.Options = options;
            elseif strcmp(obj.OptimizationConfiguration.Solver, "fminsearch")
                error(message("quantum:annealing:qaoa:SolverOptionsMustBeFminsearch"));
            else
                error(message("quantum:annealing:qaoa:SolverOptionsMustBeForSolver", ...
                    obj.OptimizationConfiguration.Solver));
            end
            
        end

    end

    methods (Hidden)
        circuit = quboAnsatz(obj, Q, angles)
        [meas, objVals] = randsampleObjective(obj, sv, Q)
    end

    methods (Hidden, Static)
        circuit = qaoaAnsatz(targetQubits, theta, numLayers, costLayerFcn, nameValuePairs)
    end

    % Custom Display methods
    methods (Access = protected)
        groups = getPropertyGroups(obj)
    end

end

function mustBeValidOptimoptionsOrStructOrEmpty(val)

if ~isempty(val) && ~(isa(val, "optim.options.SolverOptions") || isstruct(val))
    error(message("quantum:annealing:qaoa:OptimoptionsRequired"));
end

end