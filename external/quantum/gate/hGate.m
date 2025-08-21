function g = hGate(qubits)
%HGATE  Hadamard gate
%
%   g = HGATE(targetQubit) applies a Hadamard gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, HGATE returns a column
%   vector of gates, each representing an HGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                         [1 1; 1 -1] / sqrt(2).
%
%   It transforms the qubit state between the Z-basis (|0>, |1>) and the
%   X-basis (|+>, |->).
%
%   Example:
%       % Construct Hadamard gate
%       gate = hGate(1)
%       M = getMatrix(gate)
%
%       % Construct array of Hadamard gates
%       gates = hGate(1:4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, xGate, yGate, zGate,
%   idGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("h", qubits);
