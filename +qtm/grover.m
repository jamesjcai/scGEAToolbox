function [] = grover(n, x0)
%   grover  Simulate the Grover's algorithm. Print the outcome m with
%           probability p for the given input.

X = [0 1; 1 0];
H = 1 / sqrt(2) * [1 1; 1 -1];

% prepare H^n
Hn = 1;
for j = 1:n,
    Hn = kron(Hn, H);
end

% prepare a basis state bs0
bs0 = zeros(2^n, 1);
bs0(1) = 1;

% prepare the Grover operator
S0 = 2 * bs0 * bs0' - eye(2^n);
G = kron(Hn * S0 * Hn, eye(2)) * query(n, x0);

% run the Grover's algorithm
psi = G^floor(pi / 4 * 2^(n / 2)) * kron(Hn, H * X) * kron(bs0, [1 0]');

% find the most likely outcome m
p = 0;
for j = 1:2^n
    prob = psi(2 * j - 1)^2 + psi(2 * j)^2;
    if prob > p
        m = j - 1;
        p = prob;
    end
end

disp(['get outcome m = ', num2str(m), ' with probability ', num2str(p)])


end



function [O] = query(n, x0)
%   query   Return the matrix representation of the query gate for the function
%           which f(x) = 1 if x = x0, f(x) = 0 otherwise.

Of = eye(2^n);
Of(x0 + 1, x0 + 1) = -1;

O = kron(Of, eye(2));
end


% https://github.com/usami/grover/tree/master
