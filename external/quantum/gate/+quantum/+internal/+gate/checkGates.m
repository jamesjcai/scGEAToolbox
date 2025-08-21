function checkGates(gates, numQubits)
%

%   Copyright 2021-2022 The MathWorks, Inc.

if ~isa(gates, 'quantum.gate.QuantumGate')
    error(message('quantum:quantumCircuit:invalidGate'));
end

if nargin > 1
    gatesMaxQubit = quantum.internal.gate.getMaxQubit(gates);
    if gatesMaxQubit > numQubits
        error(message('quantum:quantumCircuit:invalidQubitInGate', gatesMaxQubit, numQubits));
    end
end
end
