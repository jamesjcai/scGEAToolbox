function obj = unpack(obj,idx)
%UNPACK  Unpack CompositeGate objects in a quantumCircuit.
%
%   ucirc = UNPACK(circ) returns a new quantumCircuit ucirc for which every
%   CompositeGate has been unpacked into its containing gates.
%
%   ucirc = UNPACK(circ, id) returns a new quantumCircuit ucirc where the
%   gate with index id in circ.Gates has been unpacked. circ.Gates(id)
%   must be a CompositeGate.
%
%   ucirc = UNPACK(circ, "recursive") returns a new quantumCircuit ucirc
%   for which every CompositeGate has been unpacked into its containing
%   gates recursively. The returned ucirc.Gates contains no CompositeGate
%   objects.
%
%   Example:
%       % Construct quantumCircuit and unpack CompositeGate objects.
%       cInner = quantumCircuit([cxGate(1, 2), rxGate(2, pi/3)]);
%       gates = [hGate(1)
%                compositeGate(cInner, [1 2])
%                compositeGate(cInner, [3 1])];
%       circ = quantumCircuit(gates);
%       circ.Gates
%       ucirc = unpack(circ);
%       ucirc.Gates
%       ucirc3 = unpack(circ, 3);
%       ucirc3.Gates
%
%   Example:
%       % Construct recursive quantumCircuit and unpack
%       cInner = quantumCircuit([cxGate(1, 2), rxGate(2, pi/3)]);
%       gates = [hGate(1); compositeGate(cInner, [1 2])];
%       c = quantumCircuit(gates, 3);
%       gatesOuter = [compositeGate(cInner, [3 5])
%                     hGate(3)
%                     compositeGate(c, [4 5 1])];
%       cOuter = quantumCircuit(gatesOuter);
%       cOuter.Gates
%       uc = unpack(cOuter);
%       uc.Gates
%       ucRecursive = unpack(cOuter, 'recursive');
%       ucRecursive.Gates
%
%   See also quantumCircuit, quantum.gate.CompositeGate

%   Copyright 2021-2025 The MathWorks, Inc.

if ~isscalar(obj)
    error(message('quantum:quantumCircuit:mustBeScalar'))
end

gates = obj.Gates;

if nargin == 1
    idx = 'all';
end

if matlab.internal.math.partialMatch(idx, 'all')
    % Unpack CompositeGates at the top level
    gates = replaceCompositeGates(gates);
elseif matlab.internal.math.partialMatch(idx, 'recursive')
    % Unpack all CompositeGates
    while ~isa(gates, 'quantum.gate.SimpleGate') && ~isempty(gates)
        gates = replaceCompositeGates(gates);
    end
elseif isscalar(idx) && isnumeric(idx) && isreal(idx) && floor(idx) == idx && ...
        idx >= 1 && isfinite(idx) && isa(gates(idx),'quantum.gate.CompositeGate')
    % Unpack CompositeGate at index
    childGates = extractGates(gates(idx));
    if isscalar(gates)
        % Handle the scalar case separately so output gates has the expected
        % class of children gates and does not get concatenated with the
        % empty class of the parent gates.
        gates = childGates;
    else
        gates = [gates(1:idx-1); childGates; gates(idx+1:end)];
    end
else
    error(message('quantum:quantumCircuit:unpackInvalidInput'));
end

if isempty(gates) && isa(gates, 'quantum.gate.CompositeGate')
    gates = quantum.gate.QuantumGate.empty(0,1);
end

obj.Gates = gates;

if nargout == 0
    warning(message('quantum:quantumCircuit:NotModified'));
end
end

function gates = extractGates(cg)
gates = cg.Gates;
map = cg.TargetQubits;
for ii=1:length(gates)
    qb = getQubits(gates(ii));
    qb = map(qb);
    gates(ii) = setQubits(gates(ii), qb);
end
end

function gates = replaceCompositeGates(gates)
if isscalar(gates)
    % Handle the scalar case separately so output gates has the expected
    % class of children gates and does not get concatenated with the
    % empty class of the parent gates.
    gates = extractGates(gates);
    % Gates is a QuantumGate or SimpleGate array
else
    ii = 1;
    while ii <= numel(gates)
        if ~isa(gates(ii),'quantum.gate.CompositeGate')
            ii = ii+1;
            continue;
        end
        childGates = extractGates(gates(ii));
        gates = [gates(1:ii-1); childGates; gates(ii+1:end)];
        ii = ii + numel(childGates);
    end
end
end
