function gatesOut = cancelUCREntangleGates(gatesIn, trgt, rotGate, entgGate)
% Internal use only for the ucrgates.m helper.
% The gatesIn is the local filtered gates array for the uniform-controlled
% rotation on the target qubit. It should only contain gates with types of
% the rotGate and entgGate function handles.

% Copyright 2023-2024 The MathWorks, Inc.
arguments
    gatesIn quantum.gate.SimpleGate
    trgt
    rotGate
    entgGate
end

% Indices of remaining rotation gates
rotIdx = find([gatesIn.Type]==extractBefore(func2str(rotGate),"Gate"));

gatesOut = [];
start = 1;
for ii = 1:length(rotIdx)
    % Consecutive block of two-qubit gates between rotation gates
    blk = gatesIn(start:rotIdx(ii)-1);
    if isscalar(blk)
        gatesOut = [gatesOut; blk];
    elseif ~isempty(blk)
        % Add the two-qubit gates that don't cancel out
        gatesOut = [gatesOut; entgGate(reducedCtrlIdx(blk), trgt)];
    end

    % Add the rotation gate
    gatesOut = [gatesOut; gatesIn(rotIdx(ii))];
    start = rotIdx(ii)+1;
end

% Add the two-qubit gates that don't cancel out after the last rotation gate
gatesOut = [gatesOut; entgGate(reducedCtrlIdx(gatesIn(start:end)), trgt)];

end


function cidx = reducedCtrlIdx(g)
% Return the first control indices for the qubits with an odd number of two
% qubit gates.
ctrls = [g.ControlQubits];
[numctrls,qubits] = groupcounts(ctrls.');
isOdd = mod(numctrls, 2) == 1;
cidx = qubits(isOdd);
end