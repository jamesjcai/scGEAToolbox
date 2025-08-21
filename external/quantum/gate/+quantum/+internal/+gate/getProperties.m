function [types, ctrls, trgts, angles] = getProperties(gates)
% Internal use only. Returns unpacked properties of a gates array.

% Copyright 2024 The MathWorks, Inc.
arguments
    gates quantum.gate.QuantumGate
end

if isempty(gates)
    types = string.empty;
    ctrls = double.empty(1,0);
    trgts = double.empty(1,0);
    angles = double.empty(1,0);
    return
end

if isa(gates, 'quantum.gate.SimpleGate')
    % Concatenate all gates at once
    types = [gates.Type];
    ctrls = [gates.ControlQubits];
    trgts = [gates.TargetQubits];
    angles = [gates.Angles];
else
    L = length(gates);
    types = cell(1,L);
    ctrls = cell(1,L);
    trgts = cell(1,L);
    angles = cell(1,L);
    for ii = 1:L
        g = gates(ii);
        if isa(g, 'quantum.gate.CompositeGate')
            [typesCG, ctrlsCG, trgtsCG, anglesCG] = getProperties(g);
            types{ii} = typesCG;
            ctrls{ii} = ctrlsCG;
            trgts{ii} = trgtsCG;
            angles{ii} = anglesCG;
        else
            types{ii} = g.Type;
            ctrls{ii} = g.ControlQubits;
            trgts{ii} = g.TargetQubits;
            angles{ii} = g.Angles;
        end
    end
    types = [types{:}];
    ctrls = [ctrls{:}];
    trgts = [trgts{:}];
    angles = [angles{:}];
end
end