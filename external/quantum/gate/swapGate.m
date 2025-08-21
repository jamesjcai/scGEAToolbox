function g = swapGate(targetQubit1, targetQubit2)
%SWAPGATE  Swap gate
%
%   g = SWAPGATE(targetQubit1, targetQubit2) swaps the states of the two
%   target qubits.
%
%   If target1Qubit1 or targetQubit2 is a vector of qubit indices, SWAPGATE
%   returns a column vector of gates, each representing a SWAPGATE applied
%   to the respective elements of targetQubit1 and targetQubit2.
%
%   The matrix representation of this gate applied to two target qubits
%   is
%
%                              1     0     0     0
%                              0     0     1     0
%                              0     1     0     0
%                              0     0     0     1
%
%   Example:
%       % Construct swap gate
%       gate = swapGate(1, 2)
%       M = getMatrix(gate)
%
%       % Construct array of swap gates
%       gates = swapGate(1:4, 5)
%       gates = swapGate(1:4, 8:-1:5)
%       gates = swapGate(1, 2:5)
%
%   See also quantumCircuit, quantum.gate.SimpleGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("swap", targetQubit1, targetQubit2);
