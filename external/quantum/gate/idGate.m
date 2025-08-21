function g = idGate(qubits)
%IDGATE  Identity gate
%
%   g = IDGATE(targetQubit) applies an identity gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, IDGATE returns a column
%   vector of gates, each representing an IDGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                           [1 0; 0 1].
%
%   It leaves the qubit's state unchanged.
%
%   Example:
%       % Construct identity gate
%       gate = iGate(1)
%       M = getMatrix(gate)
%
%       % Construct array of identity gates
%       gates = idGate(1:4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, xGate, yGate, zGate,
%   hGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("id", qubits);
