function g = rxGate(targetQubit, theta)
%RXGATE  X-axis rotation gate
%
%   g = RXGATE(targetQubit, theta) applies an X-axis rotation gate to a
%   target qubit. This gate rotates the qubit state around the X-axis by an
%   angle of theta.
%
%   If targetQubit or theta is a vector of qubit indices or angles, RXGATE
%   returns a column vector of gates, each representing an RXGATE applied
%   to the respective element of targetQubit using the respective angle of
%   theta.
%
%   The matrix representation of this gate applied to a qubit is
%
%              cos(theta/2)       -1i*sin(theta/2)
%              -1i*sin(theta/2)   cos(theta/2)
%
%   It performs a rotation around the X axis of the Bloch Sphere
%   representation for the target qubit state.
%
%   Example:
%       % Construct X-axis rotation gate
%       gate = rxGate(1, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of X-axis rotation gates
%       gates = rxGate(1:4, pi/3)
%       gates = rxGate(1:4, pi/4*(0:3))
%       gates = rxGate(2, pi/4*(0:3))
%
%   See also quantumCircuit, quantum.gate.SimpleGate, ryGate, rzGate,
%   r1Gate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("rx", targetQubit, theta);
