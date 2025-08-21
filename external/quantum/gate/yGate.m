function g = yGate(qubits)
%YGATE  Pauli Y gate
%
%   g = YGATE(targetQubit) applies a Pauli Y gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, YGATE returns a column
%   vector of gates, each representing a YGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                               0     -1i
%                               1i      0
%
%   Example:
%       % Construct Pauli Y gate
%       gate = yGate(1)
%       M = getMatrix(gate)
%
%       % Construct array of Pauli Y gates
%       gates = yGate(1:4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, xGate, zGate, hGate,
%   idGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("y", qubits);
