function g = tGate(qubits)
%TGATE  T gate
%
%   g = TGATE(targetQubit) applies a T gate to a target qubit.
%
%   If targetQubit is a vector of qubit indices, TGATE returns a column
%   vector of gates, each representing a TGATE applied to targetQubit(i).
%
%   The matrix representation of this gate applied to a qubit is
%
%                         [1 0; 0 (1+1i)/sqrt(2)].
%
%   It is equivalent to r1Gate(1, pi/4), and applying tGate twice to the
%   same qubit has the same effect as applying one sGate to that qubit.
%   Applying tGate four times has the same effect as applying zGate once.
%
%   See also quantumCircuit, quantum.gate.SimpleGate, tiGate, sGate, siGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("t", qubits);
