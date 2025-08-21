function A = assignmentMatrix(states, calibrations)
% Internal use only. Returns matrix A with elements aij representing the
% probability of state j being converted to state i by measurement error.

%   Copyright 2024 The MathWorks, Inc.

% References:
% [1] P. D. Nation, H. Kang, N. Sundaresan, J. M. Gambetta. "Scalable
% Mitigation of Measurement Errors on Quantum Computers." PRX Quantum,
% Nov. 2021. American Physical Society.

S = char(states)=='1';
[L, numQubits] = size(S);

A = ones(L*L, 1);
C = reshape(calibrations, 4, size(calibrations,3));
for k = 1:numQubits
    ind = S(:, k) + 2*S(:, k)' + 1;
    A = A.*C(ind, k);
end
% Columns sum to 1
A = reshape(A, L, L);
A = A./sum(A,1);
end