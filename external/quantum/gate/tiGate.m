function g = tiGate(qubits)
%TIGATE  Inverse T gate
%
%   g = TIGATE(targetQubit) applies an inverse T gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, TIGATE returns a column
%   vector of gates, each representing a TIGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                         [1 0; 0 (1-1i)/sqrt(2)].
%
%   It is equivalent to r1Gate(1, -pi/4), and applying tiGate twice to the
%   same qubit has the same effect as applying one siGate to that qubit.
%   Applying tiGate four times has the same effect as applying zGate once.
%
%   See also quantumCircuit, quantum.gate.SimpleGate, tGate, sGate, siGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("ti", qubits);
