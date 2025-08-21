function M = getMatrix(obj)
%GETMATRIX  Unitary matrix representation of a quantumCircuit.
%   M = getMatrix(circ) returns the unitary matrix representing
%   quantumCircuit circ. This matrix has size 2^N-by-2^N, where N =
%   circ.NumQubits is the number of qubits in the circuit.
%
%   Example:
%       % Construct and return the matrix representation of a quantumCircuit.
%       gates = [hGate(1); cxGate(1, 2)];
%       circ = quantumCircuit(gates);
%       M = getMatrix(circ)
%
%   See also quantumCircuit, quantumCircuit/simulate,
%   quantum.gate.SimpleGate/getMatrix, quantum.gate.CompositeGate/getMatrix

%   Copyright 2021-2022 The MathWorks, Inc.

if isinf(2^obj.NumQubits) % Can overflow, give useful error for that case.
    error(message('MATLAB:pmaxsize'));
end
M = eye(2^obj.NumQubits);

for g = 1:length(obj.Gates)
    M = applyToState(obj.Gates(g), M, obj.NumQubits);
end
