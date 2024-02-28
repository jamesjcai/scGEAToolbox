Q = [0 -1 2;...
    -1 0 4;...
    2 4 0];
c = [-5 6 -4];
d = 12;
qprob = qubo(Q,c,d);
ts = tabuSearch(Display="iter",MaxStallTime=0.01);
result = solve(qprob,Algorithm=ts);

% https://arxiv.org/ftp/arxiv/papers/1811/1811.11538.pdf
Q=[-5 2 4 0; 2 -3 1 0; 4 1 -8 5; 0 0 5 -6];
qprob=qubo(Q);
r2=solve(qprob);
r.BestX
