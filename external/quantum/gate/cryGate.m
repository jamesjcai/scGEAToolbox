function g = cryGate(controlQubit, targetQubit, theta)
%CRYGATE  Controlled Y-axis rotation gate
%
%   g = CRYGATE(controlQubit, targetQubit, theta) applies a controlled
%   Y-axis rotation gate to a target qubit based on the state of a control
%   qubit. This gate rotates the qubit state around the Y-axis by an angle
%   of theta.
%
%   If controlQubit, targetQubit, or theta is a vector of qubit indices or
%   angles, CRYGATE returns a column vector of gates, each representing a
%   CRYGATE with respective values from controlQubit, targetQubit and
%   theta.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   ryGate(targetQubit, theta).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                   1     0     0               0
%                   0     1     0               0
%                   0     0     cos(theta/2)    -sin(theta/2)
%                   0     0     sin(theta/2)    cos(theta/2)
%
%   Example:
%       % Construct controlled Y-axis rotation gate
%       gate = cryGate(1, 2, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of controlled Y-axis rotation gates
%       gates = cryGate(1:4, 5, pi/3)
%       gates = cryGate(1:4, 2:5, pi/4*(0:3))
%       gates = cryGate(2, 3:6, pi/4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, crxGate, crzGate,
%   cr1Gate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("cry", controlQubit, targetQubit, theta);
