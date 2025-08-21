function cg = qftGate(targetQubits)
%QFTGATE Quantum Fourier Transform
%
%   cg = QFTGATE(targetQubits) returns a CompositeGate cg that applies the
%   Quantum Fourier Transform to the target qubits.
%
%   See also quantumCircuit, quantum.gate.CompositeGate, compositeGate,
%   mcxGate

%  Copyright 2022 The MathWorks, Inc.

%  References:
%  [1] M. Nielsen and I. Chuang, "Quantum Computation and Quantum Information."
%  Cambridge Series on Information and the Natural Sciences, 217-219, 2010.

arguments
    targetQubits {mustBeVector(targetQubits, 'allow-all-empties'), mustBeInteger, mustBePositive}
end

if length(unique(targetQubits))~=length(targetQubits)
    error(message("quantum:gates:matchingQubits"))
end

N = length(targetQubits);
% apply the QFT algorithm
gates = [];
for target = 1:N %#ok<*AGROW>
    gates = [gates; hGate(target); cr1Gate(target+1:N, target, pi./2.^(1:(N-target)))];
end

% reverse the order of the returned qubits
f = floor(N/2);
gates = [gates; swapGate(1:f, N:-1:N-f+1)];

cg = compositeGate(quantumCircuit(gates), [targetQubits(:)]);
cg.Name = "qft";
end
