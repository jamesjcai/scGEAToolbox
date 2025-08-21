function cg = ucrxGate(controlQubits, targetQubit, angles, NameValueArgs)
%UCRXGATE Uniform controlled X-axis rotation gate
%
%   cg = UCRXGATE(controlQubits, targetQubit, angles) returns a CompositeGate
%   that applies an X-axis rotation to a target qubit for each computational
%   basis state of the control qubits. The rotation angles must be a vector
%   with length 2^L, where L is the number of control qubits, containing an
%   angle for each of the control states.
%
%   cg = UCRXGATE(___, RotationThreshold=THRES) removes rotation gates
%   with angle magnitude below the threshold THRES. The threshold must be a
%   positive real number or one of the following:
%
%   "auto" - (default) The threshold is 2*pi*eps to remove angles that are
%            approximately 0.
%   "none" - No gates are removed.
%
%   When the control qubits precede the target qubit, the matrix representation
%   is block-diagonal, where rxMi represents the 2x2 X-axis rotation matrix
%   for the input angles(i).
%
%                         [rxM1        0  ]
%                         [     .         ]
%                         [        .      ]
%                         [  0        rxMn]
%
%   See also compositeGate, unitaryGate, initGate, ucrzGate, ucryGate

%   Copyright 2023-2024 The MathWorks, Inc.

% References:
% [1] M. Möttönen, J. J. Vartiainen, V. Bergholm, and M. M. Salomaa.
% "Quantum circuits for general multiqubit gates". Phys Rev Lett. September 2004
% [2] D. Camps and R. Van Beeumen. "FABLE: Fast Approximate Quantum Circuits
% for Block-Encodings". IEEE International Conference on Quantum Computing
% and Engineering. April 2022.
arguments
    controlQubits {mustBeVector(controlQubits, 'allow-all-empties'), mustBeInteger, mustBePositive}
    targetQubit (1,1) {mustBeInteger, mustBePositive}
    angles {mustBeVector(angles, 'allow-all-empties'), mustBeReal, mustBeFinite}
    NameValueArgs.RotationThreshold {quantum.internal.gate.mustBeRotationThreshold(NameValueArgs.RotationThreshold)} = "auto"
end

qubits = [controlQubits(:); targetQubit(:)];
if length(unique(qubits))~=length(qubits)
    error(message("quantum:gates:matchingQubits"))
end

Nc = length(controlQubits);
if length(angles) ~= 2^Nc
    error(message("quantum:CompositeGate:invalidUCRAngles", 2^Nc))
end

% Set numeric value for rotation threshold
thres = NameValueArgs.RotationThreshold;
if matlab.internal.math.partialMatch(thres, "none")
    thres = 0;
elseif matlab.internal.math.partialMatch(thres, "auto")
    thres = 2*pi*eps;
end

if Nc==0
    g = rxGate(targetQubit, angles);
    cg = quantum.internal.gate.filterRotationGates(g, thres);
else
    gates = quantum.internal.gate.ucrgates(angles, @rxGate, @czGate, thres);
    cg = compositeGate(gates, [controlQubits targetQubit], Name="ucrx");
end
end