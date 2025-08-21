function gates = mcxsplitvchain(ctrls, trgt, ancilla)

% Copyright 2022 The MathWorks, Inc.

if length(ctrls)==3
    gates = quantum.internal.gate.mcxfullvchain(ctrls, trgt, ancilla);
else
    ctrls1 = ctrls(1:2:end);
    ctrls2 = ctrls(2:2:end);
    Nc1 = length(ctrls1);
    Nc2 = length(ctrls2);
    
    g1 = quantum.internal.gate.mcxfullvchain(ctrls1, ancilla, ctrls2(1:Nc1-2));
    g2 = quantum.internal.gate.mcxfullvchain([ctrls2 ancilla], trgt, ctrls1(1:Nc2-1));
    
    gates = [g1; g2; g1; g2];
end
end