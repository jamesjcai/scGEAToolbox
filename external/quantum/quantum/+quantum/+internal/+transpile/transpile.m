function [qasm, finalMap] = transpile(circuit, adj, basisMode)
%TRANSPILE  Internal use only.
%
%   [qasm, finalMap] = TRANSPILE(circuit, adj, basisMode) returns OpenQASM
%   code qasm for the circuit that satsfies device adjacency matrix adj and
%   basis gate mode basisMode. The output map finalMap assigns each circuit
%   qubit to a device qubit based on a heuristic for adding SWAP gates [1].
%
%   Adjacency element a(i,j) is true if the two-qubit basis gate of the device
%   can be applied from qubit i to qubit j. Rows and columns represent
%   control and target qubits by convention.
%   Basis mode is an integer representing the basis gate set of the device.
%   All returned code uses gates from this basis set. There are 2 modes:
%
%   5 - {ecr, rz, sx, x}
%   6 - {cz, rz, sx, x}
%
% Copyright 2024 The MathWorks, Inc.

% References:
% [1] Fei Hua, Meng Wang, Gushu Li, Bo Peng, Chenxu Liu, Muqing Zheng,
% Samuel Stein, Yufei Ding, Eddy Z. Zhang, Travis S. Humble, Ang Li.
% "QASMTrans: A QASM based Quantum Transpiler Framework for NISQ Devices."
% arXiv preprint arXiv:2308.07581 (2023)
arguments
    circuit (1,1) quantumCircuit
    adj
    basisMode (1,1) double
end

numCircuitQubits = circuit.NumQubits;

% Check adjacency and remove self loops
G = digraph(adj, 'omitselfloops');
if numCircuitQubits > numnodes(G)
    error(message("quantum:transpile:InvalidAdjacencySize"))
end
bins = conncomp(G, 'Type', 'weak');
if ~all(bins==1)
    % Each device qubit must have a path to all other qubits
    error(message("quantum:transpile:InvalidAdjacencyComponent"))
end
adj = logical(full(adjacency(G)));

[types, ctrls, trgts, angles] = getProperties(circuit);

initialMap = 0:numCircuitQubits-1;

if ~isempty(types)

    randStr = RandStream('dsfmt19937','Seed',0);
    initialMapGuess = randperm(randStr, numCircuitQubits);

    [instructions, finalMap] = quantum.internal.transpile.transpileMex( ...
        types, ctrls, trgts, angles, numCircuitQubits, adj, basisMode, initialMapGuess);
else
    instructions = "";
    finalMap = initialMap;
end

% Build the complete code
header = "OPENQASM 3.0;"+newline+...
    'include "stdgates.inc";'+newline;

creg = newline+sprintf("bit[%g] c;", numCircuitQubits)+newline+newline;

measurements = join("c["+initialMap+"] = measure $"+finalMap+";", newline);

qasm = header+...
    creg+...
    instructions+...
    measurements;
end