function g = rxxGate(targetQubit1, targetQubit2, theta)
%RXXGATE  Ising XX coupling gate
%
%   g = RXXGATE(targetQubit1, targetQubit2, theta) applies an Ising XX
%   coupling gate to two target qubits with a phase parameter theta.
%
%   If one of the target qubits or theta is a vector of qubit indices or
%   angles, RXXGATE returns a column vector of gates, each representing an
%   RXXGATE applied to the respective element of targetQubit1 and
%   targetQubit2 using the respective angle of theta.
%
%   The matrix representation of this gate applied to two qubits is
%
%           expm(-1i*theta/2*kron(X, X))   where   X = [0 1; 1 0]
%
%   X is the matrix representation of the Pauli X gate.
%
%   Example:
%       % Construct RXX gate
%       gate = rxxGate(1, 2, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of RXX gates
%       gates = rxxGate(1:4, 5, pi/3)
%       gates = rxxGate(1, 2:3, pi/4*(0:1))
%       gates = rxxGate(2:3, 4:5, pi/4*(0:1))
%
%   See also quantumCircuit, quantum.gate.SimpleGate, ryyGate, rzzGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("rxx", targetQubit1, targetQubit2, theta);
