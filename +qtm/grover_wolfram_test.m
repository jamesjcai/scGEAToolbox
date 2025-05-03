n = 3;
gates1 = [hGate(1:n-1);
         mcxGate(1:n-1, n, []);
         hGate(1:n-1);
         xGate(1:n-1);
         czGate(2, 1);
         xGate(1:n-1);
         hGate(1:n-1)];
cg = quantumCircuit(gates1);
plot(cg)
s = simulate(cg);
[states,P] = querystates(s,[1 2])

% ----------------------------

gates1 = [hGate(1:n-1);
         czGate(2, 1);
         hGate(1:n-1);
         xGate(1:n-1);
         czGate(2, 1);
         xGate(1:n-1);
         hGate(1:n-1)];
cg = quantumCircuit(gates1);
plot(cg)
s = simulate(cg);
[states,P] = querystates(s,[1 2])

%% ---------------------------------------

N = 2;
cg = diffuser1(N);
figure; plot(cg); title(cg.Name);
y = checkmatrix(cg, N)

cg = diffuser2(N);
figure; plot(cg); title(cg.Name);
y = checkmatrix(cg, N)

cg = diffuser3(N);
figure; plot(cg); title(cg.Name);
y = checkmatrix(cg, N)

cg = diffuser4(N);
figure; plot(cg); title(cg.Name);
%y = checkmatrix(cg, N)

function cg = diffuser1(n)
    gates = [hGate(1:n);
             xGate(1:n);
             hGate(n);
             mcxGate(1:n-1, n, []);
             hGate(n);
             xGate(1:n);
             hGate(1:n)];
    cg = quantumCircuit(gates,Name="Diffuser1");
end

function cg = diffuser2(n)
    gates = [hGate(1:n);
             xGate(1:n);
             czGate(1,2);
             xGate(1:n);
             hGate(1:n)];
    cg = quantumCircuit(gates,Name="Diffuser2");
end

function cg = diffuser3(n)
    gates = [hGate(1:n);
             xGate(1:n);
             hGate(2);
             cxGate(1,2);
             hGate(2);             
             xGate(1:n);
             hGate(1:n)];
    cg = quantumCircuit(gates,Name="Diffuser3");
end


function cg = diffuser4(n)
gates = [hGate(1:n);
         xGate(1:n);
         mcxGate(1:n, n+1, []);
         xGate(1:n);
         hGate(1:n)];
    cg = quantumCircuit(gates,Name="Diffuser4");
end

function y = checkmatrix(cg, N)
    e = ones(2^N,1);
    e = e/norm(e);
    M = eye(2^N) - 2*(e*e');
    y = norm(getMatrix(cg) - M);
end