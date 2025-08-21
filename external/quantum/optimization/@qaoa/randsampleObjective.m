function [meas, objVal] = randsampleObjective(obj, Q, theta)
%RANDSAMPLEOBJECTIVE Randomly sample the state and calculate objective
%
%   [MEAS, OBJVAL] = RANDSAMPLEOBJECTIVE(OBJ, Q, THETA) creates the ansatz
%   from the quadratic matrix Q and the current angles, THETA. The ansatz
%   is then simulated and the QuantumState is randomly sampled with
%   OBJ.NumShots samples. The resulting QuantumMeasurement object, MEAS and
%   the objective value at each measured state, OBJVAL, is returned.

%   Copyright 2024 The MathWorks, Inc.

% Create and simulate the QUBO Circuit
circuit = quboAnsatz(obj, Q, theta);
sv = simulate(circuit);

% Measure circuit numShots times
meas = randsample(sv, obj.NumShots);

% Calculate objective at each measured state
numStates = numel(meas.MeasuredStates);
objVal = zeros(numStates, 1);
for i = 1:numStates
    x = double(char(meas.MeasuredStates(i))=='1').';
    objVal(i) = x'*Q*x; 
end
