Q = [0 -1 2;...
    -1 0 4;...
    2 4 0];
c = [-5 6 -4];
d = 12;
qprob = qubo(Q,c,d);
ts = tabuSearch(Display="iter",MaxStallTime=0.01);
result = solve(qprob,Algorithm=ts)