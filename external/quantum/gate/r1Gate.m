function g = r1Gate(targetQubit, theta)
%R1GATE  Z-axis rotation gate with global phase
%
%   g = R1GATE(targetQubit, theta) applies a Z-axis rotation gate with
%   global phase to a target qubit. This gate shifts the phase of the
%   |1> state of the target qubit.
%
%   If targetQubit or theta is a vector of qubit indices or angles, R1GATE
%   returns a column vector of gates, each representing an R1GATE applied
%   to the respective element of targetQubit using the respective angle of
%   theta.
%
%   The matrix representation of this gate applied to a qubit is
%
%                        1        0
%                        0        exp(1i*theta)
%
%   The r1Gate is equal to the rzGate up to a global phase factor - in
%   Bloch sphere representation, both gates apply rotations around the Z
%   axis for the target qubit state.
%
%   Note that r1Gate(qubit, pi) is the same as zGate(qubit), while
%   rzGate(qubit, pi) has a different global phase. The r1Gate or rzGate
%   applied to a single qubit will result in a global phase difference that
%   is not measurable on a quantum device. But, applying the controlled
%   version of these gates, such as cr1Gate or crzGate, there is a
%   difference in relative phase that is measurable on a quantum device.
%
%   Example:
%       % Construct Z-axis rotation gate (with global phase)
%       gate = r1Gate(1, pi/5)
%       M = getMatrix(gate)
%
%       % Construct array of Z-axis rotation gates (with global phases)
%       gates = r1Gate(1:4, pi/3)
%       gates = r1Gate(1:4, pi/4*(0:3))
%       gates = r1Gate(2, pi/4*(0:3))
%
%   See also quantumCircuit, quantum.gate.SimpleGate, rxGate, ryGate,
%   rzGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.SimpleGate("r1", targetQubit, theta);
