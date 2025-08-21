function obs = qubo2ising(qp)
%QUBO2ISING Convert qubo problem to observable
%
%   obs = qubo2ising(qb) converts the qubo problem or matrix to an Ising 
%   observable. The returned observable represents the qubo objective as a
%   weighted sum of Pauli strings. Each binary variable is transformed using
%
%                          x = (I-Z)/2 
% 
%   where I and Z are the 2x2 Identity and Pauli Z matrix. The observable
%   matrix contains all qubo objective values along the diagonal.
%
%   Example:
%       % Construct from qubo object
%       qp = qubo([0 -2; -2 0], [5,1], 7);
%       obs = qubo2ising(qp);
%       % Diagonal of the observable matrix contains all objective values
%       H = getMatrix(obs);
%
%   Example:
%       % Construct from qubo matrix
%       Q = [1 3 5; 3 5 7; 5 7 9];
%       obs = qubo2ising(Q);
%
%   See also observe, observable

%   Copyright 2025 The MathWorks, Inc.  

% References:
% MathWorks. (n.d.). What is a QUBO? from 
% https://www.mathworks.com/help/matlab/math/what-is-a-qubo.html
arguments
    qp (1,1) qubo
end

% Q is real symmetric and may have diagonal elements 
Q = qp.QuadraticTerm;
[qii, qjj, qwts] = find(triu(Q, 1)/2);

% Add diagonal and quadratic contributions
c = qp.LinearTerm;
[cii, ~, cwts] = find(-(c+sum(Q,2))/2);

numQ = length(qii);
numC = length(cii);

numPauli = numQ+numC+1;
pauli = repmat('I', [numPauli qp.NumVariables]);
wts = zeros(numPauli, 1);

% Quadratic
wts(1:numQ) = qwts;
for kk = 1:numQ
    qbs = [qii(kk) qjj(kk)];
    pauli(kk, qbs) = 'Z';
end

% Linear
wts(numQ+1:end-1) = cwts;
for kk = 1:numC
    pauli(numQ+kk, cii(kk)) = 'Z';
end

% Constant represents the all-Identity Pauli string
wts(end) = sum(qwts) + sum(c+diag(Q))/2 + qp.ConstantTerm;

obs = observable(pauli, wts);
end
