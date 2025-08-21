function g = ccxGate(controlQubit1, controlQubit2, targetQubit)
%CCXGATE  Controlled-Controlled Pauli X gate
%
%   g = CCXGATE(controlQubit1, controlQubit2, targetQubit) applies a
%   controlled Pauli X gate to a target qubit based on the states of two
%   control qubits. This gate is also known as the CCNOT or Toffoli gate.
%
%   If one of the control qubits or the target qubit is a vector of qubit
%   indices, CCXGATE returns a column vector of gates, each representing a
%   CCXGATE applied to the respective elements of controlQubit1,
%   controlQubit2 and targetQubit.
%
%   If one or both of controlQubit1 and controlQubit2 are in state |0>,
%   this gate does nothing. If both controlQubit1 and controlQubit2 are in
%   state |1>, this gate is equivalent to xGate(targetQubit).
%
%   The matrix representation of this gate applied to controlQubits 1 and
%   2 and targetQubit 3 is
%
%                 1     0     0     0     0     0     0     0
%                 0     1     0     0     0     0     0     0
%                 0     0     1     0     0     0     0     0
%                 0     0     0     1     0     0     0     0
%                 0     0     0     0     1     0     0     0
%                 0     0     0     0     0     1     0     0
%                 0     0     0     0     0     0     0     1
%                 0     0     0     0     0     0     1     0
%
%   Example:
%       % Construct Controlled-Controlled X gate
%       gate = ccxGate(1, 2, 3)
%       M = getMatrix(gate)
%
%       % Construct array of Controlled-Controlled X gates
%       gates = ccxGate(1:4, 5, 6)
%       gates = ccxGate(1:4, 2:5, 3:6)
%       gates = ccxGate(1, 2, 3:5)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, cxGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("ccx", controlQubit1, controlQubit2, targetQubit);
