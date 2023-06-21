
layer1 = [ryGate(1,0);
          ryGate(2,0);
          ryGate(3,0);];

layer2 = [cryGate(1,2,0);
          cryGate(2,1,0);
          cryGate(2,3,0);
          cryGate(3,2,0);
          cryGate(1,3,0);
          cryGate(3,1,0)];

C = quantumCircuit([layer1; layer2]);
plot(C)
S = simulate(C);
S.BasisStates
S.Amplitudes
% f = formula(S)
histogram(S)
[states,P] = querystates(S)
p = probability(S,2,"1")

M = randsample(S,50)
T = table(M.Counts,M.MeasuredStates,VariableNames=["Counts","States"])



