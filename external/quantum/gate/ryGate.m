function g = ryGate(targetQubit, theta)
%RYGATE  Y-axis rotation gate
%
%   g = RYGATE(targetQubit, theta) applies a Y-axis rotation gate to a
%   target qubit. This gate rotates the qubit state around the Y-axis by an
%   angle of theta.
%
%   If targetQubit or theta is a vector of qubit indices or angles, RYGATE
%   returns a column vector of gates, each representing an RYGATE applied
%   to the respective element of targetQubit using the respective angle of
%   theta.
%
%   The matrix representation of this gate applied to a qubit is
%
%              cos(theta/2)   -sin(theta/2)
%              sin(theta/2)   cos(theta/2)
%
%   It performs a rotation around the Y axis of the Bloch Sphere
%   representation for the target qubit state.
%
%   Example:
%       % Construct Y-axis rotation gate
%       gate = ryGate(1, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of Y-axis rotation gates
%       gates = ryGate(1:4, pi/3)
%       gates = ryGate(1:4, pi/4*(0:3))
%       gates = ryGate(2, pi/4*(0:3))
%
%   See also quantumCircuit, quantum.gate.SimpleGate, rxGate, rzGate,
%   r1Gate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("ry", targetQubit, theta);
