function [gates, wasFiltered] = filterRotationGates(gates, thres)
% Internal use only.
% Filter the input gates array to remove single-qubit rotations with an
% angle less than the threshold and return if any gates were filtered. 
% All gates must be a SimpleGate since this does not filter children of a 
% CompositeGate.

% Copyright 2023-2024 The MathWorks, Inc.
arguments
    gates quantum.gate.SimpleGate
    thres
end

angles = {gates.Angles};
keepGate = false(size(angles));
for ii = 1:length(angles)
    if isempty(angles{ii}) || abs(angles{ii}) >= thres
        % Keep gates with an empty angle (non-rotation gates) and ones
        % with angles above the threshold
        keepGate(ii) = true;
    end
end

wasFiltered = ~all(keepGate);

if isscalar(keepGate) && ~keepGate
    % All gates are filtered. The scalar case is handled explicitly since
    % indexing in scalar with false returns 0x0 instead of 0x1 SimpleGate
    gates = quantum.gate.SimpleGate.empty(0,1);
else
    gates = gates(keepGate);
end
end