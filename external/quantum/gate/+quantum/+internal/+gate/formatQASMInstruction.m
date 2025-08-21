function instr = formatQASMInstruction(obj, name)

%   Copyright 2021-2023 The MathWorks, Inc.

if nargin==1
    name = obj.type;
end

if isa(obj, 'quantum.gate.CompositeGate')
    qargs = join("q["+string([obj.TargetQubits]-1)+"]", ",");
    instr = sprintf("%s %s;", name, qargs) + newline;
else
    % Expect obj to be internal gate class inheriting from quantum.internal.gate.gates.Gate in this branch.
    if ~isempty(getControlQubits(obj))
        qargs = join("q["+string([getControlQubits(obj) getTargetQubits(obj)]-1)+"]", ",");
    else
        qargs = join("q["+string(getTargetQubits(obj)-1)+"]", ",");
    end

    if ~isempty(obj.angles)
        formatstr = join(repmat("%.14g", 1, length(obj.angles)), ",");
        params = sprintf(formatstr, obj.angles);
        instr = sprintf("%s(%s) %s;", name, params, qargs) + newline;
    else
        instr = sprintf("%s %s;", name, qargs) + newline;
    end
end
