
% https://www.mathworks.com/matlabcentral/fileexchange/4711-gsa-m
%This script simulates the Quantum Mechanical Lov Grover's 
%Search Alghorithm.
clear all;

%-----parameters-----------
nqubits = 3;%number of q-bits
n = 2 ^ nqubits;%nnumber of elements in database
findmode = mod(round(n * rand + 1), n);%desired element

%-----defining quantum gates
d = - eye(n) + 2 / n;  %diffusion transform

%{
e = ones(2^N,1);
e = e/norm(e);
M = eye(2^N) - 2*(e*e');     % reflection matrix

d = - eye(2^N) + 2 / (2^N);  % diffusion transform


bs0 = zeros(2^n, 1); bs0(1) = 1;
S0 = 2 * bs0 * bs0' - eye(2^n);
Hn * S0 * Hn


%}



oracle = eye(n);  %oracle
oracle(findmode, findmode) = - 1; 


%{
    oracle = eye(8);  %oracle
    oracle(3, 3) = - 1;

    % U_w is oracle
    w = [0 0 1 0 0 0 0 0]';
    U_w = eye(8) - 2*w*w'
    % https://link.springer.com/article/10.1038/s41598-020-80153-z?fromPaywallRec=false

    isequal(U_w, oracle)
%}

%--calculate the optimal number of iterations---
finish = round(pi / 4 * sqrt( n ));

%--step(i)--initialization----
psistart = ones(n,1) / sqrt(n);
psi = psistart * exp(1i * rand);
psi0 = psi;

%step (ii)--algorithm body----
for steps = 1:finish
    %steps
    psi = d * oracle * psi;
    %probability(steps) = psi(findmode) * conj(psi(findmode));
end
%see the probability dynamics
%plot(probability);
%see the result distribution
figure;
stem(psi.*conj(psi));