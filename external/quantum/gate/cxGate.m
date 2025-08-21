function g = cxGate(controlQubit, targetQubit)
%CXGATE  Controlled Pauli X gate
%
%   g = CXGATE(controlQubit, targetQubit) applies a controlled Pauli X gate
%   to a target qubit based on the state of a control qubit. This gate is
%   equivalent to the controlled NOT gate cnotGate(controlQubit,
%   targetQubit).
%
%   If controlQubit or targetQubit is a vector of qubit indices, CXGATE
%   returns a column vector of gates, each representing a CXGATE applied to
%   the respective elements of controlQubit and targetQubit.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   xGate(targetQubit).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                              1     0     0     0
%                              0     1     0     0
%                              0     0     0     1
%                              0     0     1     0
%
%   Example:
%       % Construct Controlled X gate
%       gate = cxGate(1, 2)
%       M = getMatrix(gate)
%
%       % Construct array of Controlled X gates
%       gates = cxGate(1:4, 5)
%       gates = cxGate(1:4, 2:5)
%       gates = cxGate(1, 2:5)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, cnotGate, cyGate,
%   czGate, chGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("cx", controlQubit, targetQubit);
