function cp = plot(obj, varargin)
%PLOT  Plots a quantumCircuit
%
%   PLOT(circ) plots the quantumCircuit circ.
%
%   PLOT(___, QubitBlocks=blocks) distinguishes blocks of qubits by
%   separating them with a red dashed line in the plot. Blocks should be a
%   vector, where each element is the size of a block. The block sizes must
%   sum up to the number of qubits in circ.
%
%   PLOT(___, NumRows=nr) specifies the number of rows to use when wrapping
%   a circuit over multiple rows. By default, NumRows is determined based
%   on the input circuit.
%
%   PLOT(___, QubitLabelLocation=loc) specifies the location of labels for
%   the qubit lines as one of the following values:
%        'left'   - Qubit lines are labeled on the left.
%        'right'  - Qubit lines are labeled on the right.
%        'none'   - Qubit lines are not labeled.
%        'both'   - Qubit lines are labeled on the left and right.
%   By default, this value is chosen based on the input circuit.
%
%   PLOT(parent,circ,___) plots quantumCircuit circ in the figure, panel, or
%   tab specified by parent.
%
%   C = PLOT(...) returns a QuantumCircuitChart object. Use the methods and
%   properties of this object to inspect and adjust the plotted circuit.
%
%   Example:
%       % Construct and plot a quantumCircuit.
%       gates = [hGate(1); cxGate(1, 2)];
%       circ = quantumCircuit(gates);
%       plot(circ)
%
%   Example:
%       % Construct and plot a quantumCircuit, and show the effect of
%       % setting the QubitBlocks option:
%       gates = ccxGate([1 3 4 5], [2 6 7 8], [6 7 8 9]);
%       circ = quantumCircuit(gates);
%       plot(circ, QubitBlocks=[5 3 1])
%
%   See also quantumCircuit, quantum.gate.CompositeGate/plot

%   Copyright 2021-2025 The MathWorks, Inc.

args = varargin;
nameOffset = 1;

% Check if the first input argument is a graphics object to use as parent.
if isa(obj,'matlab.graphics.Graphics')
    % plot(parent,circ,___)
    args = [args(:),{'Parent',obj}];
    obj = args{1};
    args(1) = [];
    nameOffset = 2;
end

if ~isa(obj, 'quantumCircuit')
    error(message('quantum:quantumCircuit:plotMustBeQuantumCircuit', nameOffset));
end
if ~isscalar(obj)
    error(message('quantum:quantumCircuit:mustBeScalar'))
end

quantum.internal.gate.checkQubitBlocksForPlot(obj.NumQubits, args);

p = quantum.gate.QuantumCircuitChart('NumQubits', obj.NumQubits, 'Gates', obj.Gates, args{:});

if nargout > 0
    cp = p;
end
end