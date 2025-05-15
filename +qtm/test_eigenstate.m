g = [ryGate(1, pi/4); hGate(1); xGate(1); hGate(1)];
c = quantumCircuit(g);
plot(c)
s = simulate(c);
[states,P] = querystates(s);
formula(s)