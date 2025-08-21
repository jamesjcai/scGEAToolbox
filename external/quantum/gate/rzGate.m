function g = rzGate(targetQubit, theta)
%RZGATE  Z-axis rotation gate
%
%   g = RZGATE(targetQubit, theta) applies a Z-axis rotation gate to a
%   target qubit. This gate rotates the qubit state around the Z-axis by an
%   angle of theta.
%
%   If targetQubit or theta is a vector of qubit indices or angles, RZGATE
%   returns a column vector of gates, each representing an RZGATE applied
%   to the respective element of targetQubit using the respective angle of
%   theta.
%
%   The matrix representation of this gate applied to a qubit is
%
%              exp(-1i*theta/2)   0
%              0                  exp(1i*theta/2)
%
%   It performs a rotation around the Z axis of the Bloch Sphere
%   representation for the target qubit state.
%
%   Example:
%       % Construct Z-axis rotation gate
%       gate = rzGate(1, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of Z-axis rotation gates
%       gates = rzGate(1:4, pi/3)
%       gates = rzGate(1:4, pi/4*(0:3))
%       gates = rzGate(2, pi/4*(0:3))
%
%   See also quantumCircuit, quantum.gate.SimpleGate, r1Gate, rxGate,
%   ryGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("rz", targetQubit, theta);
