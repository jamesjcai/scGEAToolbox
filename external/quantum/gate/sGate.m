function g = sGate(qubits)
%SGATE  S gate
%
%   g = SGATE(targetQubit) applies an S gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, SGATE returns a column
%   vector of gates, each representing an SGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                         [1 0; 0 1i].
%
%   It is equivalent to r1Gate(1, pi/2), and applying sGate twice to the
%   same qubit has the same effect as applying one zGate to that qubit.
%
%   See also quantumCircuit, quantum.gate.SimpleGate, siGate, tGate, tiGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("s", qubits);
