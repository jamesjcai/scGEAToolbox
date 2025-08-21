function g = compositeGate(c, targetQubits, varargin)
%COMPOSITEGATE  Construct a composite gate for quantum computing
%
%   cg = COMPOSITEGATE(circ, targetQubits) constructs a CompositeGate from
%   quantumCircuit circ. Each qubit of circ is mapped to a qubit of an
%   outer circuit containing cg through the targetQubits vector of qubit
%   indices.
%   The length of targetQubits must match the number of qubits of circ. The
%   Name property of circ is copied over to Name property of cg.
%   
%   cg = COMPOSITEGATE(gates, targetQubits) constructs a CompositeGate from
%   an array of gates. Each qubit that the gates act on is mapped to a
%   qubit of an outer circuit containing cg through the targetQubits vector
%   of qubit indices.
%   This syntax errors if any gate uses a qubit with index larger than the
%   length of targetQubits.
%
%   cg = COMPOSITEGATE(___, Name=name) additionally sets the Name property
%   of CompositeGate cg to name.
%
%   Example:
%       % Construct a CompositeGate from a quantumCircuit:
%       innerGates = [hGate(1); cxGate(1, 2)];
%       innerCirc = quantumCircuit(innerGates, Name="bell");
%       outerGates = [hGate(1:4)
%                     compositeGate(innerCirc, [1 3])
%                     compositeGate(innerCirc, [2 4])];
%       outerCirc = quantumCircuit(outerGates);
%       plot(outerCirc)
%       % Click on a CompositeGate block in the plot - a new figure
%       % showing the internal gates of the block will appear.
%
%   See also quantumCircuit, quantum.gate.CompositeGate, mcxGate, qftGate

%   Copyright 2022 The MathWorks, Inc.

g = quantum.gate.CompositeGate(c, targetQubits, varargin{:});
