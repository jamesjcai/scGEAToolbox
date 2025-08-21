function g = siGate(qubits)
%SIGATE  Inverse S gate
%
%   g = SIGATE(targetQubit) applies an inverse S gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, SIGATE returns a column
%   vector of gates, each representing an SIGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                         [1 0; 0 -1i].
%
%   It is equivalent to r1Gate(1, -pi/2), and applying siGate twice to the
%   same qubit has the same effect as applying one zGate to that qubit.
%
%   See also quantumCircuit, quantum.gate.SimpleGate, sGate, tGate, tiGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("si", qubits);
