classdef (Sealed) QuantumState < matlab.mixin.Scalar
    %QUANTUMSTATE  State of a quantum circuit.
    %
    %   The quantumCircuit/simulate method returns a QuantumState object.
    %   Additionally, a QuantumState object can be constructed directly
    %   using:
    %
    %   s = quantum.gate.QuantumState(basisString) takes a string
    %   containing the letters "0", "1", "+" or "-" and translates it
    %   into a quantum state. Each letter represents the state of one
    %   qubit, with the first letter mapping to the first qubit.
    %
    %   s = quantum.gate.QuantumState(vec) takes a numeric vector of length
    %   2^numQubits, and normalizes it in euclidean norm to produce a valid
    %   quantum state. The Amplitudes property is set to vec / norm(vec).
    %   vec must be finite and must not be all-zero.
    %
    %   Example:
    %       % QuantumState returned from simulating a quantumCircuit
    %       gates = [hGate(1); cxGate(1, 2)];
    %       circ = quantumCircuit(gates);
    %       s = simulate(circ)
    %
    %   Example:
    %       % Directly construct a QuantumState from a string
    %       s = quantum.gate.QuantumState("01+-")
    %       formula(s)
    %       % Directly construct a QuantumState from a complex vector
    %       s = quantum.gate.QuantumState([1; -1; 1i; -1i])
    %       s.Amplitudes
    %
    %   QuantumState properties:
    %      BasisStates      - String array of basis states.
    %      Amplitudes       - Complex vector of amplitudes.
    %      NumQubits        - Number of qubits.
    %
    %   QuantumState methods:
    %      formula           - Formula representation of the state.
    %      histogram         - Histogram plot of possible states.
    %      querystates       - Query all possible states of a set of qubits.
    %      probability       - Probability of making given measurement.
    %      randsample        - Randomly sample the state.
    %
    %   See also quantum.gate.QuantumMeasurement, quantumCircuit/simulate

    %   Copyright 2021-2023 The MathWorks, Inc.

    properties(Dependent)
        %BASISSTATES - String array of basis states
        %   BasisStates is a vector of strings representing all
        %   combinations of the basis states for NumQubits qubits. Each
        %   element of the vector is a string of letters "0" and "1" with
        %   length NumQubits, and the vector length is 2^NumQubits.
        %
        %   The elements in BasisStates map to the Amplitudes property,
        %   which represents the quantum state as a linear combination of
        %   these basis states.
        %
        %   See also quantum.gate.QuantumState
        BasisStates
    end

    properties(GetAccess = public, SetAccess = private)
        %   AMPLITUDES is a complex vector of length 2^NumQubits, where
        %   each element is the amplitude of the corresponding basis state.
        %   This mapping represents the quantum state as a linear
        %   combination of the BasisStates.
        %
        %   See also quantum.gate.QuantumState
        Amplitudes

        %NUMQUBITS - Number of qubits.
        %   NUMQUBITS is the number of qubits.
        %
        %   Both BasisStates and Amplitudes properties are vectors of
        %   length 2^NumQubits.
        %
        %   See also quantum.gate.QuantumState
        NumQubits
    end

    methods
        function s = get.BasisStates(sv)
            s = statesInBasis(sv.NumQubits);
        end
    end

    methods
        % constructor
        function sv = QuantumState(state)
            if (isstring(state) && isscalar(state)) || (ischar(state) && (isrow(state) || isequal(state, '')))
                % specify state via a bit string
                state = char(state);

                if strlength(state) == 0
                    sv.Amplitudes = 1;
                    sv.NumQubits = 0;
                    return
                end
                numQubits = length(state);

                c = cell(1, numQubits);
                for ii=1:numQubits
                    if state(ii) == '0'
                        c{ii} = [1; 0];
                    elseif state(ii) == '1'
                        c{ii} = [0; 1];
                    elseif state(ii) == '+'
                        c{ii} = [1; 1]/sqrt(2);
                    elseif state(ii) == '-'
                        c{ii} = [1; -1]/sqrt(2);
                    else
                        error(message('quantum:QuantumState:invalidBitstringConstructor'));
                    end
                end

                sv.Amplitudes = quantum.internal.gate.krond(c{:});
                sv.NumQubits = numQubits;

            elseif isfloat(state) && isvector(state)
                % pass in state vector directly
                numQubits = log2(length(state));
                if numQubits ~= floor(numQubits) || numQubits < 0
                    error(message('quantum:QuantumState:invalidLength'))
                end

                state = state ./ norm(state);
                if ~all(isfinite(state))
                    error(message('quantum:QuantumState:nonfiniteVector'))
                end

                sv.Amplitudes = state(:);
                sv.NumQubits = numQubits;
            else
                error(message('quantum:QuantumState:invalidQuantumState'));
            end
        end

        function [str, states, amps] = formula(sv, NameValueArgs)
            %FORMULA Formula representation of the state
            %
            %   str = FORMULA(s) returns a string representation of the
            %   QuantumState object as a formula.
            %
            %   [str, states, amplitudes] = FORMULA(s) additionally returns
            %   a string vector of states and a numeric vector of
            %   amplitudes that are used to form the output formula in str.
            %
            %   [...] = FORMULA(__, Basis=basis) specifies the basis in
            %   which each qubit should be represented. The basis input can
            %   be "Z" ("0" and "1" in the output states), "X" ("+" and "-"
            %   in the output), a string specifying "X" or "Z" for each
            %   selected qubit, or "auto" (default).
            %   With "auto", the "Z" basis is preferred, but "X" basis is
            %   chosen for any qubit which can always be represented as "+"
            %   or "-".
            %
            %   [...] = FORMULA(s, Threshold=thr) specifies a probability
            %   threshold where a state with probability less than thr is
            %   not included in the output.
            %   Default: s.NumQubits*eps('like', s.Amplitudes)
            %   If Threshold is set to "none", output contains every possible
            %   state, even if it has probability 0 of being measured.
            %
            %   Example:
            %       % Display the result of a quantumCircuit simulation as formulas.
            %       c = quantumCircuit(cxGate(1, 2));
            %       for inState=["00", "01", "10", "11"]
            %           s = simulate(c, inState);
            %           disp("|" + inState + "> -> " + formula(s));
            %       end
            %
            %       % Use outputs from formula to programmatically verify behavior
            %       c = quantumCircuit(cxGate(1, 2));
            %       for inState=["00", "01", "10", "11"]
            %           s = simulate(c, inState);
            %           [str, outState, amps] = formula(s);
            %           assert(isequal(amps, 1));
            %           sIn = char(inState);
            %           sOut = char(outState);
            %           if sIn(1) == '0'
            %               assert(isequal(sIn, sOut));
            %           else
            %               assert(sOut(1) == sIn(1) && sOut(2) ~= sIn(2));
            %           end
            %       end
            %
            %   See also quantum.gate.QuantumState/querystates

            arguments
                sv {mustBeA(sv, 'quantum.gate.QuantumState')}
                NameValueArgs.Basis {mustBeBaseTypeAllowAuto} = 'auto'
                NameValueArgs.Threshold {quantum.internal.gate.mustBeThreshold(NameValueArgs.Threshold)} = defaultThreshold(sv)
            end

            if matlab.internal.math.partialMatch(NameValueArgs.Threshold, "none")
                NameValueArgs.Threshold = 0;
            end

            basis = NameValueArgs.Basis;

            if matlab.internal.math.partialMatch(basis, 'auto')
                basis = detectBasisTypes(sv.Amplitudes, sv.NumQubits, NameValueArgs.Threshold);
            else
                if strlength(basis) ~= 1 && strlength(basis) ~= sv.NumQubits
                    error(message("quantum:QuantumState:invalidBasisAuto"));
                end

                basis = standardizeBaseType(basis, sv.NumQubits);
            end

            states = statesInBasis(sv.NumQubits, basis);
            amps = transformAmplitudes(sv.Amplitudes, sv.NumQubits, basis);

            nzStates = abs(amps).^2 >= NameValueArgs.Threshold;
            states = states(nzStates);
            amps = amps(nzStates);

            if isreal(amps)
                ampString = string(amps);
            else
                ampString = "(" + string(amps) + ")";
            end

            str = strjoin(pad(ampString) + " * |" + states + ">", " +" + newline);
        end

        function [states, probs] = querystates(sv, qubits, NameValueArgs)
            %QUERYSTATES Query all possible states of a set of qubits.
            %
            %   [states, probabilities] = QUERYSTATES(s) gives all possible
            %   states and their probabilities.
            %
            %   [states, probabilities] = QUERYSTATES(s, qubits) gives all
            %   possible states of the selected qubits and their
            %   probabilities.
            %
            %   [states, probabilities] = QUERYSTATES(___, Threshold=thr)
            %   provides a threshold on the probabilities, below which a
            %   state is not included in the output.
            %   Default: s.NumQubits*eps('like', s.Amplitudes)
            %   If Threshold is set to "none", states contains every possible
            %   state, even if it has probability 0 of being measured.
            %
            %   [states, probabilities] = QUERYSTATES(___, Basis=basis)
            %   specifies the basis in which each qubit should be
            %   represented. The basis input can be "Z" (default, "0" and
            %   "1" in the output states), "X" ("+" and "-" in the output),
            %   or a string specifying "X" or "Z" for each selected qubit.
            %
            %   Use the histogram method (with no output argument) to plot
            %   the states and their probabilities as a histogram. For
            %   QuantumState objects, the histogram method has the same
            %   syntaxes as querystates.
            %
            %   See also quantum.gate.QuantumState/histogram,
            %   quantum.gate.QuantumMeasurement/querystates.
            arguments
                sv {mustBeA(sv, 'quantum.gate.QuantumState')}
                qubits {mustBeQubits(qubits, sv)} = 1:sv.NumQubits
                NameValueArgs.Basis {mustBeBaseType} = 'Z'
                NameValueArgs.Threshold {quantum.internal.gate.mustBeThreshold(NameValueArgs.Threshold)} = defaultThreshold(sv)
            end

            if matlab.internal.math.partialMatch(NameValueArgs.Threshold, "none")
                NameValueArgs.Threshold = 0;
            end

            qubitsUnset = isequal(qubits, 1:sv.NumQubits);

            % Process Qubits NVP
            % Process Basis NVP and apply basis transformation
            basis = NameValueArgs.Basis;

            if strlength(basis) ~= 1 && length(qubits) ~= strlength(basis)
                error(message('quantum:QuantumState:invalidBasis'))
            end

            basis = standardizeBaseType(basis, length(qubits));
            allBasis = repmat('Z', 1, sv.NumQubits);
            allBasis(qubits) = basis;

            amplitudes = transformAmplitudes(sv.Amplitudes, sv.NumQubits, allBasis);

            probs = abs(amplitudes).^2;

            % Reduce number of qubits to the requested qubits set
            if ~qubitsUnset
                probs = reduceQubits(probs, sv.NumQubits, qubits);
            end

            states = statesInBasis(length(qubits), basis);

            % Find states with probability at and above the threshold
            nzStates = probs >= NameValueArgs.Threshold;
            states = states(nzStates);
            probs = probs(nzStates);
        end

        function h = histogram(sv, varargin)
            %HISTOGRAM Histogram plot of possible states.
            %
            %   HISTOGRAM(s) plots a histogram showing each possible state
            %   and the probability of measuring the state.
            %
            %   HISTOGRAM(s, qubits) plots all possible states of the
            %   selected qubits and their probabilities.
            %
            %   HISTOGRAM(___, Threshold=thr) provides a threshold on the
            %   probabilities, below which a state is not included in the
            %   output.
            %   Default: s.NumQubits*eps('like', s.Amplitudes)
            %   If Threshold is set to "none", every possible state is plotted,
            %   even if it has probability 0 of being measured.
            %
            %   HISTOGRAM(___, Basis=basis) specifies the basis in which
            %   each qubit should be represented. The basis input can be
            %   "Z" (default, "0" and "1" in the plotted output), "X" ("+"
            %   and "-" in the output), or a string specifying "X" or "Z"
            %   for each selected qubit.
            %
            %   h = HISTOGRAM(...) also returns a Histogram object. Use
            %   this to inspect and adjust properties of the histogram.
            %
            %   Use querystates (with two output arguments) to return the
            %   possible states and their probabilities. For QuantumState
            %   objects, the querystates method has the same syntaxes as
            %   histogram.
            %
            %   See also quantum.gate.QuantumState/querystates,
            %   quantum.gate.QuantumMeasurement/histogram.

            % Same NVPs as querystates
            [k, p] = querystates(sv, varargin{:});

            hInner = histogram('Categories', categorical(k), 'BinCounts', p);
            ylabel('Probability')
            if nargout > 0
                h = hInner;
            end
        end

        function p = probability(sv, qubits, state)
            %PROBABILITY Probability of making given measurement
            %
            %   p = probability(s, qubits) returns the probability of
            %   having all specified qubits being measured in state "1".
            %   This probability is based on the probability distribution
            %   of the quantumState s. Output p is a scalar between 0 and 1.
            %
            %   p = PROBABILITY(s, qubits, state) additionally specifies
            %   which state these qubits should be measured as. Input state
            %   can be "0", "1", "+", or "-" meaning all qubits are being
            %   measured in one of these states. Input state can also be a
            %   string with length(qubits) of such characters to specify
            %   the measured state of each qubit. The default state is "1".
            %
            %   See also quantum.gate.QuantumMeasurement/probability.

            arguments
                sv {mustBeA(sv, 'quantum.gate.QuantumState')}
                qubits {mustBeQubits(qubits, sv)}
                state {mustBeBasisString} = '1'
            end

            state = char(state);
            if ~isscalar(state) && length(state) ~= length(qubits)
                error(message("quantum:QuantumState:invalidBitstring"));
            end

            if isscalar(state)
                state = repmat(state, 1, length(qubits));
            end

            basis = repmat('Z', size(state));
            basis(ismember(state, '+-')) = 'X';
            state01 = state;
            state01(state=='+') = '0';
            state01(state=='-') = '1';

            basis = standardizeBaseType(basis, length(qubits));
            allBasis = repmat('Z', 1, sv.NumQubits);
            allBasis(qubits) = basis;

            amplitudes = transformAmplitudes(sv.Amplitudes, sv.NumQubits, allBasis);

            probs = abs(amplitudes).^2;

            probs = reduceQubits(probs, sv.NumQubits, qubits);

            ind = bin2dec(state01) + 1;

            p = probs(ind);
        end

        function counts = randsample(sv, numShots)
            %RANDSAMPLE Randomly sample the state
            %
            %   m = RANDSAMPLE(s, numShots) randomly samples the
            %   QuantumState s with numShots samples, and returns the
            %   resulting QuantumMeasurement object.
            %
            %   This function uses the current random stream, similary to
            %   the function rand. This means that a different result is
            %   returned on each call. Reset the random number stream to
            %   control this behavior.
            %
            %   See also quantum.gate.QuantumMeasurement, rand, rng

            arguments
                sv {mustBeA(sv, 'quantum.gate.QuantumState')}
                numShots {mustBeScalarOrEmpty, mustBeInteger} = 100
            end

            cumProbabilities = [0; cumsum(abs(sv.Amplitudes).^2)];

            samples = rand(numShots, 1);

            ind = discretize(samples, cumProbabilities);

            [indUnique, ~, bin] = unique(ind);

            count = accumarray(bin, 1, [numel(indUnique), 1]);

            allStates = statesInBasis(sv.NumQubits);
            states = allStates(indUnique);

            if isempty(count)
                % NumQubits information is lost if states is an empty
                % string array
                states = char.empty(0, sv.NumQubits);
            end

            counts = quantum.gate.QuantumMeasurement(states, count);
        end
    end

    methods (Access=private)
        function thresh = defaultThreshold(sv)
            % Multiply by NumQubits assuming that the number of gates
            % applied to a state is roughly proportional to the number of
            % qubits.
            thresh = max(1, sv.NumQubits)*eps('like', sv.Amplitudes);
        end
    end

    properties(Constant, Access='private')
        % Version of the serialization and deserialization
        % format. This is used for managing forward compatibility. Value is
        % saved in 'versionSavedFrom' when an instance is serialized.
        %
        %   1.0 : original shipping version
        version = 1.0;
    end
    methods (Hidden)
        function s = saveobj(qs)
            s = struct('Amplitudes', qs.Amplitudes, ...
                'versionSavedFrom', quantum.gate.QuantumState.version, ...
                'minCompatibleVersion', 1);
        end
    end
    methods(Hidden, Static)
        function qs = loadobj(s)
            if quantum.gate.QuantumState.version < s.minCompatibleVersion
                id = 'quantum:QuantumState:IncompatibleVersion';
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.QuantumState', getString(message(id))));
                warning(id, loadWarningString);

                qs = quantum.gate.QuantumState(1);
                return
            end

            try
                qs = quantum.gate.QuantumState(s.Amplitudes);
            catch err
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.QuantumState', err.message));
                warning(err.identifier, loadWarningString);

                qs = quantum.gate.QuantumState(1);
                return
            end
        end
    end
end

function probs = reduceQubits(probs, numQubits, qubits)
sz = 2*ones(1, numQubits);
sz(end+1:2) = 1;
probs = reshape(probs, sz);

% Reshape probs to 2x2x...x2 tensor
newsz = 2*ones(1, numQubits);
newsz(end+1:2) = 1;
probs = reshape(probs, newsz);

% Sum over all dimensions not matching one of the qubits
removeQubits = setdiff(1:numQubits, qubits);
keepDim = flip(numQubits+1-qubits);
removeDim = numQubits+1-removeQubits;
probs = sum(probs, [removeDim numQubits+1]); % last one to avoid empty dim edge-case in SUM

% Reshape back to a vector, now of length 2^length(qubits)
permvec = [keepDim removeDim];
if length(permvec) > 1
    probs = permute(probs, permvec);
end
probs = probs(:);
end


function k = statesInBasis(numQubits, basis)
% Example basis: 'XZX', meaning the first and last qubit is in
% the '+'/'-' basis, and the middle qubit is in the '0'/'1'
% basis.
if numQubits == 0
    k = "";
else
    k = dec2bin(0:2^numQubits-1, numQubits);
    if nargin > 1 && ~all(basis == 'Z')
        for jj=1:numQubits
            if basis(jj) == 'X'
                for ii=1:size(k, 1)
                    if k(ii, jj) == '0'
                        k(ii, jj) = '+';
                    else
                        k(ii, jj) = '-';
                    end
                end
            end
        end
    end
    k = string(k);
end
end

function state = transformAmplitudes(state, numQubits, basis)
% Example basis: 'XZX', meaning the first and last qubit is in
% the '+'/'-' basis, and the middle qubit is in the '0'/'1'
% basis.
H = [1 1; 1 -1]/sqrt(2);
for ii=1:numQubits
    if basis(ii) == 'X'
        state = quantum.internal.gate.applyMatToDim(state, H, numQubits-ii+1);
    end
end
end

function basis = detectBasisTypes(state, numQubits, tol)
% Transform from 'Z' basis to 'X' basis
H = [1 1; 1 -1]/sqrt(2);
for ii=1:numQubits
    state = quantum.internal.gate.applyMatToDim(state, H, numQubits-ii+1);
end
state = reshape(state, [2*ones(1, numQubits) 1 1]);
% Detect for each dimension whether it's all-zero in either '+' or '-'
basis = repmat('Z', 1, numQubits);
state = abs(state).^2;
for ii=1:numQubits
    dim = numQubits-ii+1;
    v = sum(state, [1:dim-1 dim+1:numQubits+1]); % numQubits+1 to work around empty dim behavior
    if any(v < tol)
        basis(ii) = 'X';
    end
end

end

function mustBeBaseType(str)
% Passing QuantumState instead of just numQubits as an FAV workaround
mustBeTextScalar(str)
str = lower(char(str));
if ~all(ismember(str, 'xz'))
    error(message("quantum:QuantumState:invalidBasis"));
end
end

function mustBeBaseTypeAllowAuto(str)
% Passing QuantumState instead of just numQubits as an FAV workaround
mustBeTextScalar(str)
if matlab.internal.math.partialMatch(str, 'auto')
    return
end
str = lower(char(str));
if ~all(ismember(str, 'xz'))
    error(message("quantum:QuantumState:invalidBasisAuto"));
end
end

function str = standardizeBaseType(str, numQubits)
str = char(upper(str));
if strlength(str) == 1
    str = repmat(str, 1, numQubits);
end
end

function mustBeBasisString(str)
% Passing QuantumState instead of just numQubits as an FAV workaround
mustBeTextScalar(str)
str = char(str);
if ~all(ismember(str, '01+-'))
    error(message("quantum:QuantumState:invalidBitstring"));
end
end

function mustBeQubits(qubits, sv)
numQubits = sv.NumQubits;
if ~isnumeric(qubits) || ~isvector(qubits) || ~isreal(qubits) || any(floor(qubits) ~= qubits) || ...
        numel(unique(qubits)) ~= numel(qubits) || any(qubits < 1) || any(qubits > numQubits)
    error(message("quantum:QuantumState:invalidQubits"));
end
end
