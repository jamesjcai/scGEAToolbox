function g = crzGate(controlQubit, targetQubit, theta)
%CRZGATE  Controlled Z-axis rotation gate
%
%   g = CRZGATE(controlQubit, targetQubit, theta) applies a controlled
%   Z-axis rotation gate to a target qubit based on the state of a control
%   qubit. This gate rotates the qubit state around the Z-axis by an angle
%   of theta.
%
%   If controlQubit, targetQubit, or theta is a vector of qubit indices or
%   angles, CRZGATE returns a column vector of gates, each representing a
%   CRZGATE with respective values from controlQubit, targetQubit and
%   theta.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   rzGate(targetQubit, theta).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                   1     0     0                  0
%                   0     1     0                  0
%                   0     0     exp(-1i*theta/2)   0
%                   0     0     0                  exp(1i*theta/2)
%
%   Example:
%       % Construct controlled Z-axis rotation gate
%       gate = crzGate(1, 2, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of controlled Z-axis rotation gates
%       gates = crzGate(1:4, 5, pi/3)
%       gates = crzGate(1:4, 2:5, pi/4*(0:3))
%       gates = crzGate(2, 3:6, pi/4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, crxGate, cryGate,
%   cr1Gate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("crz", controlQubit, targetQubit, theta);
