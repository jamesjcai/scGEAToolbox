function cg = initGate(targetQubits, state, NameValueArgs)
%INITGATE Gate to initialize state
%
%   cg = INITGATE(targetQubits, state) returns a CompositeGate that
%   initializes the target qubits to the input state up to a global phase.
%   The target qubits must all be in the state |0> when initGate is applied
%   for the input state to be initialized.
%
%   Input state must be a QuantumState object, basis string, or numeric
%   vector representing the same number of qubits as the number of target
%   qubits. The basis string and numeric vector must be valid inputs to
%   the QuantumState constructor. When input state is a numeric vector,
%   it is normalized as state/norm(state).
%
%   cg = INITGATE(___, RotationThreshold=THRES) removes rotation gates
%   with angle magnitude below the threshold THRES. The threshold must be a
%   positive real number or one of the following:
%
%   "auto" - (default) The threshold is 2*pi*eps to remove angles that are
%            approximately 0.
%   "none" - No gates are removed.
%
%   The CompositeGate uses O(2^N) gates to initialize the N target qubits,
%   or a small O(N) gates when input state is a basis string. It is not
%   guaranteed to use the minimal number of gates possible for a given
%   input state.
%
%   See also quantum.gate.QuantumState, ucrxGate, ucryGate, ucrzGate

%   Copyright 2023-2024 The MathWorks, Inc.

% References:
% [1] V. V. Shende, S. S. Bullock and I. L. Markov. "Synthesis of quantum-
% logic circuits". IEEE Transactions on Computer-Aided Design of Integrated
% Circuits and Systems. June 2006
arguments
    targetQubits {mustBeVector(targetQubits, 'allow-all-empties'), mustBeInteger, mustBePositive}
    state
    NameValueArgs.RotationThreshold {quantum.internal.gate.mustBeRotationThreshold(NameValueArgs.RotationThreshold)} = "auto"
end

numQubits = length(targetQubits);
if length(unique(targetQubits)) ~= numQubits
    error(message("quantum:gates:matchingQubits"))
end

if isa(state, 'quantum.gate.QuantumState')
    stateObj = state;
else
    % Validate input state using constructor
    stateObj = quantum.gate.QuantumState(state);
end

if numQubits~=stateObj.NumQubits
    error(message("quantum:CompositeGate:initIncorrectNumQubits", numQubits, stateObj.NumQubits))
end

if isempty(targetQubits)
    cg = compositeGate(quantum.gate.SimpleGate.empty(0,1), targetQubits, Name="init");
    return
end

if isstring(state) || ischar(state)
    % String state is initialized using a few simple gates
    circuit = initBasisString(state);
else
    % Vector state is initialized using uniform-controlled rotations to
    % entangle each qubit one-by-one.

    % Set numeric value for rotation threshold
    thres = NameValueArgs.RotationThreshold;
    if matlab.internal.math.partialMatch(thres, "none")
        thres = 0;
    elseif matlab.internal.math.partialMatch(thres, "auto")
        thres = 2*pi*eps;
    end

    state = stateObj.Amplitudes;
    circuit = initVector(state, thres);
end

cg = compositeGate(circuit, targetQubits, Name="init");

end

%% Helper Functions

function c = initBasisString(state)
state = char(state);
numQubits = length(state);
if all(state=='0')
    c = quantumCircuit(numQubits);
    return
end

gates = [xGate(strfind(state, ("1"|"-")))
    hGate(strfind(state,("+"|"-")))];

c = quantumCircuit(gates, numQubits);
end


function c = initVector(state, thres)

gates = [];

mag = abs(state);
ang = angle(state);
numQubits = log2(length(state));
for trgt = numQubits:-1:1
    % Compute angles to move the target qubit from where it is to |0>.
    % Each iteration cuts the magnitude and angle length in half by
    % disentangling the qubits one-by-one.
    [mag, ang, thetas, phis] = disentangle(mag, ang);

    % The ucryGate vector is equivalent when reversed and used to possibly
    % cancel gates
    rucry = flip(quantum.internal.gate.ucrgates(thetas, @ryGate, @cxGate, thres));

    ctrls = 1:trgt-1;
    if norm(phis,1) > 0
        % The input state is complex and also requires a ucrzGate
        ucrz = quantum.internal.gate.ucrgates(phis, @rzGate, @cxGate, thres);
        if ~isempty(ucrz) && ~isempty(rucry) && ...
                (ucrz(end).Type=="cx" && isequal(ucrz(end),rucry(1)))
            % Cancel out two adjacent cxGates
            ucrz(end) = [];
            rucry(1) = [];
        end
        disentg = [compositeGate(ucrz, [ctrls trgt], Name="ucrz")
            compositeGate(rucry, [ctrls trgt], Name="ucry")];
    else
        % The input state is real and only a ucryGate is required.
        disentg = compositeGate(rucry, [ctrls trgt], Name="ucry");
    end

    % Add the block of gates to disentangle the target qubit
    gates = [gates; disentg];
end

% Inverse the disentangle gates to initialize the input state
c = inv(quantumCircuit(gates, numQubits));
end

function [magOut, angleOut, nthetas, nphis] = disentangle(magIn, angleIn)
% Compute angles for ucryGate and ucrzGate to disentangle the bottom
% qubit from the input state vector of N total qubits (Theorem 9 [1]).

% The magnitude and angle components of the state vector are reshaped
% into the 2^(N-1) components of the N-1 control states of the uniform
% rotation

% ucryGate angles
magIn = reshape(magIn, 2, []);
nthetas = -2*atan2(magIn(2,:), magIn(1,:));
magOut = vecnorm(magIn, 2, 1);

% ucrzGate angles
angleIn = reshape(angleIn, 2, []);
nphis = (angleIn(1,:) - angleIn(2,:));
angleOut = (angleIn(1,:) + angleIn(2,:))/2;
end