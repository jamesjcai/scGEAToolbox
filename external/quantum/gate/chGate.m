function g = chGate(controlQubit, targetQubit)
%CHGATE  Controlled Hadamard gate
%
%   g = CHGATE(controlQubit, targetQubit) applies a controlled Hadamard
%   gate to a target qubit based on the state of a control qubit.
%
%   If controlQubit or targetQubit is a vector of qubit indices, CHGATE
%   returns a column vector of gates, each representing a CHGATE applied to
%   the respective elements of controlQubit and targetQubit.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   hGate(targetQubit).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                     1     0     0           0
%                     0     1     0           0
%                     0     0     1/sqrt(2)   1/sqrt(2)
%                     0     0     1/sqrt(2)   -1/sqrt(2)
%
%   Example:
%       % Construct Controlled Hadamard gate
%       gate = chGate(1, 2)
%       M = getMatrix(gate)
%
%       % Construct array of Controlled Hadamard gates
%       gates = chGate(1:4, 5)
%       gates = chGate(1:4, 2:5)
%       gates = chGate(1, 2:5)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, cxGate, cyGate,
%   czGate, cnotGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("ch", controlQubit, targetQubit);
