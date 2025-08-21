function g = cyGate(controlQubit, targetQubit)
%CYGATE  Controlled Pauli Y gate
%
%   g = CYGATE(controlQubit, targetQubit) applies a controlled Pauli Y gate
%   to a target qubit based on the state of a control qubit.
%
%   If controlQubit or targetQubit is a vector of qubit indices, CYGATE
%   returns a column vector of gates, each representing a CYGATE applied to
%   the respective elements of controlQubit and targetQubit.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   yGate(targetQubit).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                              1     0     0     0
%                              0     1     0     0
%                              0     0     0     -1i
%                              0     0     1i    0
%
%   Example:
%       % Construct Controlled Y gate
%       gate = cyGate(1, 2)
%       M = getMatrix(gate)
%
%       % Construct array of Controlled Y gates
%       gates = cyGate(1:4, 5)
%       gates = cyGate(1:4, 2:5)
%       gates = cyGate(1, 2:5)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, cxGate, czGate,
%   cnotGate, chGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("cy", controlQubit, targetQubit);
