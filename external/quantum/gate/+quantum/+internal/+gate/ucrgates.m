function gates = ucrgates(angles, rotGate, entgGate, thres)
% Internal helper used by ucrxGate, ucryGate, and ucrzGate that returns the
% local gates array for implementing the uniform-controlled rotation

% Copyright 2023-2024 The MathWorks, Inc.
arguments
    angles
    rotGate function_handle
    entgGate function_handle
    thres
end

% Construct local uniform-controlled rotation gate circuit with the control
% qubits at the top and target at the bottom.
Nc = log2(length(angles));
trgt = Nc+1;

if Nc==0
    gates = quantum.internal.gate.filterRotationGates(rotGate(trgt, angles), thres);
    return
end

% Apply discrete Walsh-Hadamard transformation
angles = quantum.internal.gate.dwht(angles);

% Permute the angles according to the Gray code.
code = quantum.internal.gate.buildGrayCode(Nc);
angles = angles(bin2dec(code)+1);

% Control qubits are the indices where each row of the Gray code
% differs. The first and last rows differ at the first qubit.
difIdx = code(1:2^Nc-1,:)~=code(2:2^Nc,:);
[ctrls, ~] = find(difIdx.');
ctrls(end+1) = 1;

% Merge the rotation gates between the two-qubit gates
gates = reshape([rotGate(trgt, angles)'; entgGate(ctrls, trgt)'], [], 1);

% Filter rotation gates
[gates, wasFiltered] = quantum.internal.gate.filterRotationGates(gates, thres);

if wasFiltered
    % Remove leftover two-qubit gates that cancel out.
    gates = quantum.internal.gate.cancelUCREntangleGates(gates, trgt, rotGate, entgGate);
end
end