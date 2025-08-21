function [bestx, bestf, exitflag, output] = search(obj, Q, ~)
%SEARCH Search for solutions to a QUBO
%
%   [BESTX, BESTF, EXITFLAG, OUTPUT] = SEARCH(OBJ, C, CONSTANTFORDISPLAY)
%   searches for solutions to the supplied QUBO by calling a quantum
%   approximate optimization algorithm. The iterative display reflects the
%   problem constant via the addition of CONSTANTFORDISPLAY.
%
%   Note: BESTF should not include the CONSTANTFORDISPLAY. This is added
%   back in the ABSTRACTALGORITHM superclass.

%   Copyright 2024 The MathWorks, Inc.

startTime = tic;
numLayers = obj.NumLayers;
boundValue = pi;
bound = repmat(boundValue, 2*numLayers, 1);
% Handle initial angle vectors
x0Theta = obj.InitialAngles;
if isempty(x0Theta)
    x0Theta = -boundValue + 2*boundValue*rand(2, numLayers);
end
if ~isequal(size(x0Theta), [2 numLayers])
    error(message("quantum:annealing:qaoa:InitialAnglesIncorrectSize"));
end

% Call optimization solver
objFcn = @(theta)expectedObjectiveValue(obj, theta, Q);
switch obj.OptimizationSolver
    case "surrogateopt"     
        if isempty(obj.OptimizationSolverOptions)
            opts = optimoptions("surrogateopt");
        else
            opts = obj.OptimizationSolverOptions;
        end
        opts.InitialPoints = x0Theta;
        [theta, bestExpectedValue, exitflag, outputClassical] = ...
            surrogateopt(objFcn, -bound, bound, [], [], [], [], [], opts);
    case "fmincon"
        [theta, bestExpectedValue, exitflag, outputClassical] = ...
            fmincon(objFcn, x0Theta, [], [], [], [], -bound, bound, [], obj.OptimizationSolverOptions);
    case "fminsearch"
        [theta, bestExpectedValue, exitflag, outputClassical] = ...
            fminsearch(objFcn, x0Theta, obj.OptimizationSolverOptions);        
end

% Simulate the circuit with the optimized angles
circuit = quboAnsatz(obj, Q, theta);
xstates = simulate(circuit);

% Find best X and function value
[xstr, pr] = querystates(xstates);
x = double(char(xstr)=='1').';
numX = size(x, 2);
fVals = zeros(1, numX);
for i = 1:numX
    fVals(i) = x(:, i)'*Q*x(:, i);
end
[bestf, idxBest] = min(fVals);
bestx = x(:, idxBest);

% Generate algorithm results structure
output.AllX = x;
output.AllFunctionValues = fVals;
output.BestX = bestx;
output.BestX = bestf;
output.AllProbabilities = pr';
output.NumShots = obj.NumShots;
output.BestCircuitAngles = theta;
output.ClassicalOutput = outputClassical;
output.ClockTime = toc(startTime);
output.BestExpectedValue = bestExpectedValue;
output.Circuit = circuit;

end

function expVal = expectedObjectiveValue(obj, theta, Q)

% Measure circuit numShots times and calculate objective at each state
[meas, objVals] = randsampleObjective(obj, Q, theta);

% Expectation
sumObjVals = sum(meas.Counts.*objVals);
expVal = sumObjVals/obj.NumShots;

end



