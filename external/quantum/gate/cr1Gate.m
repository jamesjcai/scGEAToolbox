function g = cr1Gate(controlQubit, targetQubit, theta)
%CR1GATE  Controlled Z-axis rotation gate with global phase
%
%   g = CR1GATE(controlQubit, targetQubit, theta) applies a controlled
%   Z-axis rotation gate with global phase to a target qubit based on the
%   state of a control qubit. This gate shifts the phase of the |1> state
%   of the target qubit.
%
%   If controlQubit, targetQubit, or theta is a vector of qubit indices or
%   angles, CR1GATE returns a column vector of gates, each representing a
%   CR1GATE with respective values from controlQubit, targetQubit and
%   theta.
%
%   If controlQubit is in state |0>, this gate does nothing. If
%   controlQubit is in state |1>, this gate is equivalent to
%   r1Gate(targetQubit, theta).
%
%   The matrix representation of this gate applied to controlQubit 1 and
%   targetQubit 2 is
%
%                   1     0     0     0
%                   0     1     0     0
%                   0     0     1     0
%                   0     0     0     exp(1i*theta)
%
%   Switching the target and control qubit in a cr1Gate has no effect on
%   the behavior.
%
%   Example:
%       % Construct controlled Z-axis rotation gate (with global phase)
%       gate = cr1Gate(1, 2, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of controlled Z-axis rotation gates (with global
%       % phases)
%       gates = cr1Gate(1:4, 5, pi/3)
%       gates = cr1Gate(1:4, 2:5, pi/4*(0:3))
%       gates = cr1Gate(2, 3:6, pi/4)
%
%   See also quantumCircuit, quantum.gate.SimpleGate, crxGate, cryGate,
%   crzGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("cr1", controlQubit, targetQubit, theta);
