A = qtm.permn([0 1], 3);
B = join(string(A), "", 2);

s = quantum.gate.QuantumState('000');
% s = quantum.gate.QuantumState([1 0 0 0 0 0 0 0]);
C = s.BasisStates;

isequal(B, C)