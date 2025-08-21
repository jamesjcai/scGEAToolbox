function y = dwht(x)
% Discrete Walsh-Hadamard Transform
% This is used to transform angles for the construction of the
% uniform-controlled rotation gates. Note this helper is equivalent to
%
% y = fwht(x, length(x), 'hadamard');
%
% from the Signal Processing Toolbox. However, this implementation is less
% efficient since it constructs and applies the full Hadamard matrix.

% References:
% [1] M. Möttönen, J. J. Vartiainen, V. Bergholm, and M. M. Salomaa.
% "Quantum circuits for general multiqubit gates". Phys Rev Lett. September 2004

% Copyright 2023-2024 The MathWorks, Inc.
arguments
    x (:,1)
end
L = length(x);
% Solve the system of angles from Equation 4 [1]
Mk = hadamard(L);
y = (1/L) * (Mk.' * x);
end