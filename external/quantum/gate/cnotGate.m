function g = cnotGate(controlQubit, targetQubit)
%CNOTGATE  Controlled Pauli X gate
%
%   g = CNOTGATE(controlQubit, targetQubit) applies a Pauli X gate to qubit
%   targetQubit, controlled by the value of controlQubit. The type of the
%   resulting gate will be "cx", calling CNOTGATE and CXGATE are equivalent.
%
%   If controlQubit or targetQubit is a vector, CNOTGATE returns a column
%   vector of gates, each representing a CNOTGATE applied to the respective
%   elements of controlQubit and targetQubit.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   xGate(targetQubit).
%
%   The matrix representation of this gate (applied to controlQubit 1 and
%   targetQubit 2) is
%
%                              1     0     0     0
%                              0     1     0     0
%                              0     0     0     1
%                              0     0     1     0
%
%   Example:
%       % Construct Controlled X gate
%       gate = cnotGate(1, 2)
%       M = getMatrix(gate)
%
%       % Construct array of Pauli X gates
%       gates = cnotGate(1:4, 5)
%       gates = cnotGate(1:4, 2:5)
%       gates = cnotGate(1, 2:5)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, cxGate, cyGate,
%   czGate, chGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("cx", controlQubit, targetQubit);
