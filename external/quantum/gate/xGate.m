function g = xGate(qubits)
%XGATE  Pauli X gate
%
%   g = XGATE(targetQubit) applies a Pauli X gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, XGATE returns a column
%   vector of gates, each representing an XGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is 
%
%                               0     1
%                               1     0
%
%   It swaps the amplitudes of the |0> and |1> states for the selected
%   qubit.
%
%   Example:
%       % Construct Pauli X gate
%       gate = xGate(1)
%       M = getMatrix(gate)
%
%       % Simulate applying X gate to a quantum circuit with a specified
%       % initial state
%       g = xGate(1);
%       circ = quantumCircuit(g, 2);
%       s = quantum.gate.QuantumState([1; -1; 1i; -1i]);
%       so = simulate(circ, s);
%       table(s.BasisStates, s.Amplitudes, so.Amplitudes, ...
%             VariableNames=["BasisStates", "Input", "Output"])
%
%       % Construct array of Pauli X gates
%       gates = xGate(1:4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, yGate, zGate, hGate,
%   idGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("x", qubits);
