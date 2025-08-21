function out = observe(state, obs, NameValueArgs)
%OBSERVE Expected value of an observable
%
%   expval = OBSERVE(state) returns the expected value of measuring the
%   quantum state in the computational basis. The expected value is a scalar
%   real number representing the average measurement result. The input can 
%   also be a quantumCircuit that prepares the state assuming all qubits
%   are initially in |0>.
%
%   expval = OBSERVE(state, obs) returns the expected value of measuring the
%   input in the basis specified by the observable using
%
%                               <s|O|s>
%
%   where |s> and O are respectively the state vector and matrix 
%   representation of the observable that act on the same number of qubits.
%
%   expval = OBSERVE(..., NumShots=nshots) approximates the expected value
%   using the specified number of random samples.
%
%   expval = OBSERVE(..., ProbabilityThreshold=alpha) approximates the
%   expected value using the conditional value-at-risk that only includes
%   the lowest expectation values up to the total probability specified by
%   alpha. This can only be set when the observable has combinations of "I" 
%   and "Z" basis measurements. 
%
%   Example:
%       % Exact expectation of a qubit in the computational basis
%       up = quantumCircuit(1);
%       mid = quantumCircuit(hGate(1));
%       down = quantumCircuit(xGate(1));
%       expval = observe(up); % Returns 1
%       expval = observe(mid); % Returns 0
%       expval = observe(down); % Returns -1
%
%   Example:
%       % Exact expectation of qubits in the basis specified by observable
%       circ = quantumCircuit(xGate([1 3]));
%       obs = observable(["ZII" "IIZ"], [3 4]);
%       expval = observe(circ, obs); % Returns -7
%
%   Example:
%       % Approximate expectation of a qubit in the specified basis using
%       % 1000 samples
%       circ = quantumCircuit([rxGate(1,pi/3); ryGate(1,pi/2)]);
%       obs = observable(["X" "Y" "Z"], 1:3);
%       expval = observe(circ, obs, NumShots=1000);
%
%   Example:
%       % Approximate expectation of qubits using 25% of smallest weighted 
%       % computational basis measurements
%       circ = quantumCircuit(hGate(1:2));
%       obs = observable(["ZI" "IZ" "ZZ"], [-1 -1 -1]);
%       expval = observe(circ, obs, ProbabilityThreshold=0.25); % Returns -3
%
%   See also observable

%   Copyright 2025 The MathWorks, Inc.

% References:
% [1] Barkoutsos, Panagiotis Kl, et al. "Improving variational quantum
% optimization using CVaR." Quantum 4 (2020).
% [2] Kohda, Masaya. et al. "Quantum expectation-value estimation by
% computational basis sampling." Phys. Rev. Research (2022).
arguments
    state (1,1)
    obs (1,1) observable = defaultObservable(state)
    NameValueArgs.NumShots (1,1) {mustBeNumeric, mustBePositive} = Inf
    NameValueArgs.ProbabilityThreshold (1,1) {mustBeInRange(NameValueArgs.ProbabilityThreshold,0,1,'exclude-lower')} = 1;
end

if isa(state, 'quantumCircuit')
    state = simulate(state);
end

numQubits = state.NumQubits;
if numQubits~=obs.NumQubits
    error(message("quantum:observable:observeIncorrectNumQubits"))
end

if ~obs.IsDiagonal && NameValueArgs.ProbabilityThreshold < 1
    error(message("quantum:observable:observeInvalidProbabilityThreshold"))
end

isExact = isinf(NameValueArgs.NumShots) && NameValueArgs.ProbabilityThreshold==1;

if isExact
    % Performance heuristic based on number of qubits and Pauli strings.
    K = length(obs.Paulis);
    useFullHermitian = ( (numQubits<=8) || (2^numQubits/numQubits<=K) );
    a = state.Amplitudes;
    if useFullHermitian
        H = getMatrix(obs);
        % Ensure full output for empty case
        out = full(a'*H*a);
    else
        % Map each Pauli to its matrix
        pU = cat(3, [0 1; 1 0], [0 -1j; 1j 0], [1 0; 0 -1]);
        X = applyPaulis(obs, a, pU, 'XYZ');
        out = obs.Weights.'*squeeze(pagemtimes(a', X));
    end
    % Remove complex round-off
    out = real(out);
else
    if obs.IsDiagonal
        out = diagEstimate(obs, state, NameValueArgs);
    else
        out = nonDiagEstimate(obs, state, NameValueArgs);
    end
end
end

%% Helper Functions

function out = diagEstimate(obs, state, NameValueArgs)
% Returns approximate expectation value from diagonal observable

if isinf(NameValueArgs.NumShots)
    [bits, p] = querystates(state);
else
    meas = randsample(state, NameValueArgs.NumShots);
    [bits, p] = querystates(meas);
end

% Spin matrix has +1/-1 for 0/1 measurements
S = (-1).^(char(bits)=='1');

% Mask matrix only considers Z observables
M = char(obs.Paulis) == 'Z';

% Reshape for broadcasting
w = obs.Weights;
N = obs.NumQubits;
S3D = reshape(S, [length(p) 1 N]);
M3D = reshape(M, [1 length(w) N]);

% Elements (i,j) of P are prod(S(i,:).^M(j,:)) which scale weights by +1 or
% -1.
P = prod(S3D.^M3D, 3);
v = P*w;

out = quantum.internal.gate.cvar(p, v, NameValueArgs.ProbabilityThreshold);
end


function out = nonDiagEstimate(obs, state, NameValueArgs)
% Returns approximate expectation value from non-diagonal observable.

% Map each Pauli to the matrix that rotates it to Z
zU = cat(3, [1 1; 1 -1], [1 -1j; 1 1j]) / sqrt(2);
X = applyPaulis(obs, state.Amplitudes, zU, 'XY');

v = zeros(size(X,3),1);
paulis = obs.Paulis;
for kk = 1:length(paulis)
    pchar = char(paulis(kk));
    iXYZ = find(pchar~='I');
    x = quantum.gate.QuantumState(X(:,:,kk));
    meas = randsample(x, NameValueArgs.NumShots);
    [bits, p] = querystates(meas, iXYZ);
    % Determine eigenvalues by parity check
    eigvals = (-1).^count(bits, "1");
    v(kk) = p.'*eigvals;
end
out = obs.Weights.'*v;
end


function X = applyPaulis(obs, s, U, map)
% Return each state vector from applying the Paulis
paulis = obs.Paulis;
N = obs.NumQubits;
numTerms = length(paulis);
X = zeros(size(s,1), 1, numTerms, class(s));
for ii = 1:numTerms
    x = s;
    pchar = char(paulis(ii));
    for jj = 1:N
        kk = strfind(map, pchar(jj));
        if ~isempty(kk)
            x = quantum.internal.gate.applyMatToDim(x, U(:,:,kk), N+1-jj);
        end
    end
    X(:,:,ii) = x;
end
end

function obs = defaultObservable(state)
if ~(isa(state, 'quantumCircuit') || isa(state, 'quantum.gate.QuantumState'))
    error(message("quantum:observable:observeInvalidInput"))
end
obs = observable(repmat('Z', [1 state.NumQubits]), 1);
end