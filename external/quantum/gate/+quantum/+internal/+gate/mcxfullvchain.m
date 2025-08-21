function gates = mcxfullvchain(ctrls, trgt, ancilla)

% Copyright 2022 The MathWorks, Inc.

Nc = length(ctrls);
Na = length(ancilla);

if Nc==2
    gates = ccxGate(ctrls(1), ctrls(2), trgt);
else
    gates = [ccxGate(ctrls(1), ctrls(2), ancilla(1))
             ccxGate(ctrls(3:Na+1), ancilla(1:Na-1), ancilla(2:Na))
             ccxGate(ctrls(end), ancilla(end), trgt)
             ccxGate(ctrls(Na+1:-1:3), ancilla(Na-1:-1:1), ancilla(Na:-1:2))
             ccxGate(ctrls(1), ctrls(2), ancilla(1))
             ccxGate(ctrls(3:Na+1), ancilla(1:Na-1), ancilla(2:Na))
             ccxGate(ctrls(end), ancilla(end), trgt)
             ccxGate(ctrls(Na+1:-1:3), ancilla(Na-1:-1:1), ancilla(Na:-1:2))];
end