function g = ryyGate(targetQubit1, targetQubit2, theta)
%RYYGATE  Ising YY coupling gate
%
%   g = RYYGATE(targetQubit1, targetQubit2, theta) applies an Ising YY
%   coupling gate to two target qubits with a phase parameter theta.
%
%   If one of the target qubits or theta is a vector of qubit indices or
%   angles, RYYGATE returns a column vector of gates, each representing an
%   RYYGATE applied to the respective element of targetQubit1 and
%   targetQubit2 using the respective angle of theta.
%
%   The matrix representation of this gate applied to two qubits is
%
%         expm(-1i*theta/2*kron(Y, Y))   where   Y = [0 -1i; 1i 0]
%
%   Y is the matrix representation of the Pauli Y gate.
%
%   Example:
%       % Construct RYY gate
%       gate = ryyGate(1, 2, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of RYY gates
%       gates = ryyGate(1:4, 5, pi/3)
%       gates = ryyGate(1, 2:3, pi/4*(0:1))
%       gates = ryyGate(2:3, 4:5, pi/4*(0:1))
%
%   See also quantumCircuit, quantum.gate.SimpleGate, rxxGate, rzzGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("ryy", targetQubit1, targetQubit2, theta);
