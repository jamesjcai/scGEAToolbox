function g = rzzGate(targetQubit1, targetQubit2, theta)
%RZZGATE  Ising ZZ coupling gate
%
%   g = RZZGATE(targetQubit1, targetQubit2, theta) applies an Ising ZZ
%   coupling gate to two target qubits with a phase parameter theta.
%
%   If one of the target qubits or theta is a vector of qubit indices or
%   angles, RZZGATE returns a column vector of gates, each representing an
%   RZZGATE applied to the respective element of targetQubit1 and
%   targetQubit2 using the respective angle of theta.
%
%   The matrix representation of this gate applied to two qubits is
%
%           expm(-1i*theta/2*kron(Z, Z))   where   Z = [1 0; 0 -1]
%
%   Z is the matrix representation of the Pauli Z gate.
%
%   Example:
%       % Construct RZZ gate
%       gate = rzzGate(1, 2, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of RZZ gates
%       gates = rzzGate(1:4, 5, pi/3)
%       gates = rzzGate(1, 2:3, pi/4*(0:1))
%       gates = rzzGate(2:3, 4:5, pi/4*(0:1))
%
%   See also quantumCircuit, quantum.gate.SimpleGate, rxxGate, ryyGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("rzz", targetQubit1, targetQubit2, theta);
