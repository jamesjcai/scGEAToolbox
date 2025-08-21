function g = czGate(controlQubit, targetQubit)
%CZGATE  Controlled Pauli Z gate
%
%   g = CZGATE(controlQubit, targetQubit) applies a controlled Pauli Z gate
%   to a target qubit based on the state of a control qubit.
%
%   If controlQubit or targetQubit is a vector of qubit indices, CZGATE
%   returns a column vector of gates, each representing a CZGATE applied to
%   the respective elements of controlQubit and targetQubit.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   zGate(targetQubit).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                              1     0     0     0
%                              0     1     0     0
%                              0     0     1     0
%                              0     0     0    -1
%
%   Switching the target and control qubit in a czGate has no effect on the
%   behavior. For this reason, the controlled Pauli Z gate is represented
%   symmetrically in plots, with a dot for both the control and the target
%   qubit.
%
%   Example:
%       % Construct Controlled Z gate
%       gate = czGate(1, 2)
%       M = getMatrix(gate)
%
%       % Construct array of Controlled Z gates
%       gates = czGate(1:4, 5)
%       gates = czGate(1:4, 2:5)
%       gates = czGate(1, 2:5)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, cxGate, cnotGate,
%   cyGate, chGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("cz", controlQubit, targetQubit);
