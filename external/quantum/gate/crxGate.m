function g = crxGate(controlQubit, targetQubit, theta)
%CRXGATE  Controlled X-axis rotation gate
%
%   g = CRXGATE(controlQubit, targetQubit, theta) applies a controlled
%   X-axis rotation gate to a target qubit based on the state of a control
%   qubit. This gate rotates the qubit state around the X-axis by an angle
%   of theta.
%
%   If controlQubit, targetQubit, or theta is a vector of qubit indices or
%   angles, CRXGATE returns a column vector of gates, each representing a
%   CRXGATE with respective values from controlQubit, targetQubit and
%   theta.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   rxGate(targetQubit, theta).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                   1     0     0                  0
%                   0     1     0                  0
%                   0     0     cos(theta/2)       -1i*sin(theta/2)
%                   0     0     -1i*sin(theta/2)   cos(theta/2)
%
%   Example:
%       % Construct controlled X-axis rotation gate
%       gate = crxGate(1, 2, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of controlled X-axis rotation gates
%       gates = crxGate(1:4, 5, pi/3)
%       gates = crxGate(1:4, 2:5, pi/4*(0:3))
%       gates = crxGate(2, 3:6, pi/4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, cryGate, crzGate,
%   cr1Gate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("crx", controlQubit, targetQubit, theta);
