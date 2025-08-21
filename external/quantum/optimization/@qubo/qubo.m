classdef (Sealed) qubo
%qubo Quadratic Unconstrained Binary Optimization 
%
%   qubo is a class that stores and solves QUBO problems
%   
%   Copyright 2022-2023 The MathWorks, Inc.    
        properties
                QuadraticTerm (:,:)  {mustBeDouble, mustBeReal, mustBeFinite}
                LinearTerm (:,1)  {mustBeDouble, mustBeReal, mustBeFinite}
                ConstantTerm (1,1)  {mustBeDouble, mustBeReal, mustBeFinite}
        end
        properties (Dependent, SetAccess = private)
                NumVariables
        end
        properties (Hidden, SetAccess = private, GetAccess = public)
            QUBOVersion
        end
        methods
            function obj = qubo(Q, c, const)
                % check inputs
                arguments
                    Q (:,:) {mustBeDouble, mustBeReal, mustBeSquare, mustBeNonempty, mustBeFinite}
                    c (:,1) {mustBeDouble, mustBeComformable(Q,c), mustBeFinite} = [];
                    const (1,1) {mustBeDouble, mustBeReal, mustBeFinite} = 0;
                end
                obj.QuadraticTerm = Q;
                obj.LinearTerm = c;
                obj.ConstantTerm = const;
                obj.QUBOVersion = 1;
            end

            function solution = solve(obj, NameValueArgs)

                arguments
                    obj
                    NameValueArgs.Algorithm (1,1) {mustBeA(NameValueArgs.Algorithm, 'quantum.internal.optimization.AbstractAlgorithm')} = tabuSearch; 
                end

                alg = NameValueArgs.Algorithm; 

                aResult = alg.run(obj);

                solution = quboResult(aResult, obj);
            end

            function obj = set.QuadraticTerm(obj, Q)
                if isempty(Q) || size(Q,1) ~= size(Q,2)
                    error(message('quantum:annealing:QUBO:QuadraticTermNotSquare'));
                elseif isempty(obj.QuadraticTerm) || all(size(Q) == size(obj.QuadraticTerm))
                    obj.QuadraticTerm = makeQSymmetric(Q);
                else
                    error(message('quantum:annealing:QUBO:QuadraticTermSizeCannotChange',size(obj.QuadraticTerm,1)));
                end
            end

            function obj = set.LinearTerm(obj, c)
                if isempty(c)
                    obj.LinearTerm = zeros(obj.NumVariables, 1); %#ok
                elseif numel(c) == obj.NumVariables %#ok
                    obj.LinearTerm = c;
                else
                    error(message('quantum:annealing:QUBO:LinearTermSizeIncorrect'));
                end
            end

            function obj = set.ConstantTerm(obj, const)
                obj.ConstantTerm = const;
            end

            function value = get.NumVariables(obj)
                value = size(obj.QuadraticTerm,1);
            end

            function fval = evaluateObjective(obj, x)
                arguments
                    obj
                    x {mustBeNumeric, mustBeDouble, mustBeBinary}
                end
                if size(x, 1) ~= obj.NumVariables
                    error(message('quantum:annealing:QUBO:xMustHaveNumVariablesRows'));
                end
                if ~ismatrix(x)
                    error(message('quantum:annealing:QUBO:xMustBeAMatrix'));
                end
                fval = NaN(1, size(x,2));
                for i = 1:size(x,2)
                    fval(i) = x(:,i)'*obj.QuadraticTerm*x(:,i) + obj.LinearTerm'*x(:,i) + obj.ConstantTerm;
                end
            end
        end
end

function Q = makeQSymmetric(Q)
    if ~issymmetric(Q)
        warning(message('quantum:annealing:QUBO:QuadraticTermNotSymmetric'));
        Q = (Q + Q.')/2;
    end
end

function mustBeSquare(Q)
    numOfRows = size(Q,1);
    numOfCols = size(Q,2);
    if (numOfRows ~= numOfCols)
        throwAsCaller(MException(message('quantum:annealing:QUBO:mustBeSquareMatrix')))
    end
end

function mustBeComformable(Q, f)
    numOfCols = size(Q,2);
    if (~isempty(f) && numel(f) ~= numOfCols)
            throwAsCaller(MException(message('quantum:annealing:QUBO:mustBeConformant')))
    end
end

function mustBeDouble(input)

    if ~isa(input, "double")
        throwAsCaller(MException(message('quantum:annealing:QUBO:mustBeDouble','QUBO')));
    end

end

function mustBeBinary(input)

if any(input ~= 0 & input ~= 1)
    throwAsCaller(MException(message('quantum:annealing:QUBO:mustBeBinary')));
end

end



