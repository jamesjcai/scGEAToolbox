function g = zGate(qubits)
%ZGATE  Pauli Z gate
%
%   g = ZGATE(targetQubit) applies a Pauli Z gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, ZGATE returns a column
%   vector of gates, each representing a ZGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                               1     0
%                               0     -1
%
%   It flips the sign of the amplitude of the |1> state for the selected
%   qubit.
%
%   Example:
%       % Construct Pauli Z gate
%       gate = zGate(1)
%       M = getMatrix(gate)
%
%       % Simulate applying Z gate to a quantum circuit with a specified
%       % initial state
%       g = zGate(1);
%       circ = quantumCircuit(g, 2);
%       s = quantum.gate.QuantumState([1; -1; 1i; -1i]);
%       so = simulate(circ, s);
%       table(s.BasisStates, s.Amplitudes, so.Amplitudes, ...
%             VariableNames=["BasisStates", "Input", "Output"])
%
%       % Construct array of Pauli Z gates
%       gates = zGate(1:4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, xGate, yGate, hGate,
%   idGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("z", qubits);
