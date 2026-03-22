function cg = mcxGate(controlQubits, targetQubit, ancillaQubits)
%MCXGATE Multi-controlled X gate
%
%   cg = MCXGATE(controlQubits, targetQubit, ancillaQubits) applies a
%   multi-controlled Pauli X gate. This gate operates on a single target
%   qubit based on the states of the control qubits, with a number of
%   ancilla qubits that don't affect the operation. If all control qubits
%   are in state |1>, then this gate applies the Pauli X gate to the target
%   qubit.
%
%   In general, the returned object is a CompositeGate that consists of a
%   number of internal gates. Each qubit of the internal gates is mapped to
%   a qubit of an outer circuit containing the composite gate. For simple
%   cases (length(controlQubits) <= 2), the returned object is a SimpleGate
%   without mapping.
%
%   The number of ancillaQubits determines the implementation of the
%   internal gates to construct the returned CompositeGate:
%
%       1 ancilla qubit (recommended)
%                       Requires a linear number of internal gates
%                       w.r.t. the number of control qubits.
%
%       n-2 ancilla qubits (where n is the number of control qubits)
%                       Requires the minimal amount of internal gates
%                       among these options (the amount still grows
%                       linearly w.r.t. the number of control qubits,
%                       but at a smaller factor)
%
%       0 ancilla qubits
%                       Requires an exponential amount of internal
%                       gates w.r.t. the number of control qubits.
%
%   The ancilla qubits can be in any state when passed to mcxGate, and are
%   returned in the same state they were originally in.
%
%   See also quantumCircuit, quantum.gate.CompositeGate, compositeGate,
%   qftGate

%  Copyright 2022 The MathWorks, Inc.

%  References:
%  [1] Barenco et al., "Elementary gates for quantum computation."
%  American Physical Society. Physical Review A, 52, 17-20, 1995.

arguments
    controlQubits {mustBeVector(controlQubits, 'allow-all-empties'), mustBeInteger, mustBePositive}
    targetQubit (1,1) {mustBeInteger, mustBePositive}
    ancillaQubits {mustBeVector(ancillaQubits, 'allow-all-empties'), mustBeInteger, mustBePositive}
end

qubits = [controlQubits(:); targetQubit(:); ancillaQubits(:)];
if length(unique(qubits))~=length(qubits)
    error(message("quantum:gates:matchingQubits"))
end

Nc = length(controlQubits);
Na = length(ancillaQubits);

if Nc==0
    cg = xGate(targetQubit);
elseif Nc==1
    cg = cxGate(controlQubits(1), targetQubit);
elseif Nc==2
    cg = ccxGate(controlQubits(1), controlQubits(2), targetQubit);
else
    ctrls = 1:Nc;
    trgt = Nc+1;
    switch Na
        case 0
            g = quantum.internal.gate.mcxgraycode(ctrls, trgt);
        case 1
            g = quantum.internal.gate.mcxsplitvchain(ctrls, trgt, Nc+2);
        otherwise
            % shorten ancillaQubits to match the number used when
            % constructing the CompositeGate
            if Na<Nc-2
                g = quantum.internal.gate.mcxsplitvchain(ctrls, trgt, Nc+2);
                ancillaQubits = ancillaQubits(1);
            else
                g = quantum.internal.gate.mcxfullvchain(ctrls, trgt, Nc+2:2*Nc-1);
                ancillaQubits = ancillaQubits(1:Nc-2);
            end
    end
    qubits = [controlQubits(:); targetQubit; ancillaQubits(:)];
    cg = compositeGate(quantumCircuit(g, length(qubits)), qubits);
    cg.Name = "mcx";
end

end
