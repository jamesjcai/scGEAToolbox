function circuit = quboAnsatz(obj, Q, theta)
%QUBOANSATZ Cost-mixer circuit for QUBO in QAOA
%
%   CIRCUIT = QUBOANSATZ(OBJ, Q, PARAMS) returns a quantumCircuit used to
%   estimate the QUBO, Q. This circuit contains numVar qubits which
%   are first put in an equal superposition via Hadamard gates. The circuit
%   then contains numLayers of gates. Each layer contains two sets of
%   gates, a cost layer and a mixer layer. The cost layer approximates the
%   cost function. The mixer layer is easy to initialize in a known ground
%   state.

%   Copyright 2024 The MathWorks, Inc.

% Map to Ising form
[J, h] = qubo2ising(Q);

% Determine the control/target qubits for the quadratic terms
[problem.Ctrl, problem.Trgt] = find(J);
problem.J = J;
problem.h = h;

% Create circuit
numVar = size(Q, 1);
circuit = qaoa.qaoaAnsatz(1:numVar, theta, obj.NumLayers, ...
    @(targetQubits, gamma)costLayer(targetQubits, gamma, problem));

end

function gates = costLayer(~, gamma, problem)

% Unwrap problem structure
ctrl = problem.Ctrl;
trgt = problem.Trgt;
J = problem.J;
h = problem.h;

% Linear terms
gates = [];
for i = 1:numel(h)        
    if h(i) ~= 0        
        gates = [gates; rzGate(i, 2*h(i)*gamma)]; %#ok
    end
end

% Quadratic terms
for ei = 1:numel(ctrl)
    gates = [gates; rzzGate(ctrl(ei), trgt(ei), 2*J(ctrl(ei), trgt(ei))*gamma)]; %#ok
end

end

function [J, h, c] = qubo2ising(Q)

% Linear term
f = diag(Q);

% Quadratic term
N = size(Q, 1);
Q(1:N+1:N^2) = 0;

% Map to Ising form
J = Q/4;
h = 0.5*(sum(Q, 2) + f);
c = 0.5*sum(f) + 0.25*sum(Q(:));

% As Q is assumed symmetric, J is symmetric and we can make J upper
% triangular. This will reduce the number of RZZ gates for the quadratic
% terms.
J = 2*triu(J);

end

