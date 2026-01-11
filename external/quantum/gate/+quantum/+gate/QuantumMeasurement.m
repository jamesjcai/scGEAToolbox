classdef (Sealed) QuantumMeasurement
    %QUANTUMMEASUREMENT Measurement of a quantum circuit
    %
    %   Running a quantumCircuit remotely on a quantum device results in a
    %   QuantumMeasurement. Simulating a quantumCircuit locally, and
    %   calling randsample on the result also returns a QuantumMeasurement.
    %
    %   Example:
    %       % Simulate a circuit locally and call randsample on the result:
    %       gates = [hGate(1); cxGate(1, 2)];
    %       circ = quantumCircuit(gates);
    %       qState = simulate(circ)
    %       qMeasure = randsample(qState, 100)
    %
    %   Example:
    %       % Run a circuit remotely on a quantum device and retrieve its
    %       % result
    %       % NOTE: This requires you to have access to AWS Braket, and will
    %       % send a quantum task that will be charged to your AWS account.
    %       gates = [hGate(1); cxGate(1, 2)];
    %       circ = quantumCircuit(gates);
    %       dev = quantum.backend.QuantumDeviceAWS("Lucy");
    %       task = run(circ, dev)
    %       task.Status
    %       wait(task)
    %       measurement = fetchOutput(task)
    %
    %   QuantumMeasurement properties:
    %      MeasuredStates   - String array of states that were measured.
    %      Counts           - Count of how often each state was measured.
    %      Probabilities    - Estimated measurement probability for each state.
    %      NumQubits        - Number of qubits.
    %
    %   QuantumMeasurement methods:
    %      histogram         - Histogram plot of measured states.
    %      querystates       - Query all measured states of a set of qubits.
    %      probability       - Probability of a specified measurement result.
    %
    %   See also quantum.gate.QuantumState, quantumCircuit/run,
    %   quantum.gate.QuantumState/randsample, quantumCircuit/simulate

    %   Copyright 2022-2025 The MathWorks, Inc.

    properties(Dependent, GetAccess=public, SetAccess=private)
        %MEASUREDSTATES - String array of states that were measured.
        %   MeasuredStates is a vector of strings representing all states
        %   that were measured at least once. It has the same length as
        %   Counts, which gives the number of times each state was measured.
        %
        %   MeasuredStates is a column string vector, where each string has
        %   NumQubits letters "0" or "1" representing whether each qubit is
        %   in |0> or |1> state.
        %
        %   See also quantum.gate.QuantumMeasurement
        MeasuredStates
    end

    properties(GetAccess=public, SetAccess=private)
        %Counts - Count of how often each state was measured
        %   Counts is a vector of the same length as MeasuredStates, giving
        %   a count of how often each state in MeasuredStates was measured.
        %
        %   Elements of COUNTS are all positive integers. If no COUNTS are
        %   returned by the remote device provider, all elements are NaN instead.
        %
        %   See also quantum.gate.QuantumMeasurement
        Counts

        %Probabilities - Estimated measurement probability for each state.
        %   Probabilities is a vector of the same length as MeasuredStates,
        %   giving the estimated measurement probability for each state in MeasuredStates.
        %
        %   Elements of PROBABILITIES are all real numbers that are not
        %   guaranteed to be positive or sum to 1 (depending on the
        %   measurement provided by the remote device provider).
        %
        %   See also quantum.gate.QuantumMeasurement
        Probabilities

        %NUMQUBITS - Number of qubits.
        %   NUMQUBITS is the number of qubits.
        %
        %   See also quantum.gate.QuantumMeasurement
        NumQubits
    end

    properties(Access=private)
        MeasuredStatesChar
    end

    methods
        function s = get.MeasuredStates(obj)
            if size(obj.MeasuredStatesChar, 1) == 0
                s = string.empty(0, 1);
            else
                s = string(obj.MeasuredStatesChar);
            end
        end

        function obj = QuantumMeasurement(states, measurements, type)
            arguments
                states {mustBeStringVectorOrCharMatrix}
                measurements (:,1) {mustBeNumeric, mustBeReal}
                type {mustBeMember(type, {'probs', 'counts'})} = 'counts'
            end

            if ischar(states)
                numEntries = size(states, 1);
                numQubits = size(states, 2);
            else
                numEntries = length(states);
                if isempty(states)
                    % just define it somehow - call constructor with char matrix input to be specific.
                    numQubits = 0;
                else
                    numQubits = strlength(states(1));
                    if ~all(strlength(states(:)) == numQubits)
                        error(message('quantum:QuantumMeasurement:sameLength'))
                    end
                end
                states = char(states(:));
                states = reshape(states, numEntries, numQubits);
            end
            assert(size(states, 2) == numQubits)

            % Usually states is unique and sorted, however we don't enforce
            % that here, instead it is done by each of our functions before
            % calling the QuantumMeasurement constructor.

            if ~all(ismember(states(:), '01'))
                error(message('quantum:QuantumMeasurement:invalidStates'))
            end

            if isequal(type, 'counts')
                mustBeInteger(measurements)
                mustBePositive(measurements)
                obj.Counts = measurements;
                obj.Probabilities = measurements ./ sum(measurements);
            else
                obj.Probabilities = measurements;
                obj.Counts = NaN(size(measurements));
            end

            if length(measurements) ~= numEntries
                error(message("quantum:QuantumMeasurement:dimMismatch"))
            end

            obj.MeasuredStatesChar = states;
            obj.NumQubits = numQubits;
        end

        function [states, probs] = querystates(obj, qubits, NameValueArgs)
            %QUERYSTATES Query all measured states of a set of qubits.
            %
            %   [states, probabilities] = QUERYSTATES(m) gives all measured
            %   states and their estimated probabilities.
            %
            %   [states, probabilities] = QUERYSTATES(m, qubits) gives all
            %   measured states of the selected qubits and their
            %   estimated probabilities.
            %
            %   [states, probabilities] = QUERYSTATES(___, Threshold=thr)
            %   provides a threshold on the probabilities, below which a
            %   measured state is not included in the output.
            %   Default: m.NumQubits*eps
            %   If Threshold is set to "none", states contains every measured
            %   state, even if it was never measured.
            %
            %   The histogram method has the same syntaxes as querystates
            %   and plots the two outputs as a histogram.
            %
            %   Use the histogram method (with no output argument) to plot
            %   the measured states and their probabilities as a histogram.
            %   For QuantumMeasurement objects, the histogram method has
            %   the same syntaxes as querystates.
            %
            %   See also quantum.gate.QuantumMeasurement/histogram,
            %   quantum.gate.QuantumState/querystates.

            arguments
                obj {mustBeScalarQuantumMeasurement(obj)}
                qubits {mustBeQubits(qubits, obj)} = 1:obj.NumQubits
                NameValueArgs.Threshold {quantum.internal.gate.mustBeThreshold(NameValueArgs.Threshold)} = defaultThreshold(obj)
            end
            
            states = obj.MeasuredStatesChar;
            probs = obj.Probabilities;

            % Reduce number of qubits to the requested qubits set
            states = states(:, qubits);
            [states, ~, id] = unique(states, 'rows');
            probs = accumarray(id, probs);

            if matlab.internal.math.partialMatch(NameValueArgs.Threshold, "none")
                NameValueArgs.Threshold = -Inf;
            end

            if (NameValueArgs.Threshold <= 0) && ~isempty(qubits)
                % Fill in non-existent states (main use case is making a
                % histogram of all possible integers, not just measured
                % ones)

                allStates = dec2bin(0:2^length(qubits)-1, length(qubits));
                if isempty(probs)
                    allProbs = nan(2^length(qubits), 1);
                else
                    allProbs = zeros(2^length(qubits), 1);
                    ind = bin2dec(states)+1;
                    allProbs(ind) = probs;
                end
                states = allStates;
                probs = allProbs;
            end

            % Return thresholded states for the qubits set
            thresStates = ~(probs < NameValueArgs.Threshold);
            states = states(thresStates, :);
            probs = probs(thresStates);
            states = string(states);
        end

        function h = histogram(obj, varargin)
            %HISTOGRAM Histogram plot of measured states.
            %
            %   HISTOGRAM(m) plots a histogram showing each measured state
            %   and the estimated probability of measuring the state.
            %
            %   HISTOGRAM(m, qubits) plots all measured states of the
            %   selected qubits and their estimated probabilities.
            %
            %   HISTOGRAM(___, Threshold=thr) provides a threshold on the
            %   probabilities, below which a measured state is not included
            %   in the output.
            %   Default: m.NumQubits*eps
            %   If Threshold is set to "none", every measured state with
            %   nonnegative probability is plotted.
            %
            %   h = HISTOGRAM(...) also returns a Histogram object. Use
            %   this to inspect and adjust properties of the histogram.
            %
            %   Use querystates (with two output arguments) to return the
            %   measured states and their probabilities. For
            %   QuantumMeasurement objects, the querystates method has the
            %   same syntaxes as histogram.
            %
            %   See also quantum.gate.QuantumMeasurement/querystates,
            %   quantum.gate.QuantumState/histogram.

            if ~isscalar(obj)
                error(message("quantum:QuantumMeasurement:mustBeScalar"))
            end

            % Same NVPs as querystates
            [states, probs] = querystates(obj, varargin{:});

            if any(probs < 0)
                % Some queried states have negative probability and are not
                % included in histogram
                warning(message("quantum:QuantumMeasurement:histogramFilteredStates"))
            end

            % Only show the non-negative states in the histogram
            nnStates = probs >= 0;
            states = states(nnStates);
            probs = probs(nnStates);

            hInner = histogram('Categories', categorical(states), 'BinCounts', probs);
            ylabel('Probability')
            if nargout > 0
                h = hInner;
            end
        end

        function p = probability(obj, qubits, state)
            %PROBABILITY Probability of a specified measurement result.
            %
            %   p = probability(m, qubits) returns the probability of
            %   having all specified qubits being measured in state "1".
            %   This probability is based on a probability estimation of
            %   the result in QuantumMeasurement m. Output p is an array
            %   with the same size as m with all elements between 0 and 1.
            %
            %   p = PROBABILITY(m, qubits, state) additionally specifies
            %   which state these qubits should be measured as. Input state
            %   can be "0" or "1", meaning all qubits are measured to be in
            %   |0> or |1> state. Input state can also be a string with
            %   length(qubits) of such characters to specify the measured
            %   state of each qubit. The default state is "1".
            %
            %   See also quantum.gate.QuantumState/probability.

            arguments
                obj {mustBeA(obj,'quantum.gate.QuantumMeasurement')}
                qubits {mustBeQubits(qubits, obj)}
                state {mustBeBasisString} = '1'
            end
            
            if ~isscalar(obj)
                sz = size(obj);
                % Use cell to handle datatypes
                p = cell(sz);
                for ii = 1:numel(obj)
                    p{ii} = probability(obj(ii), qubits, state);
                end
                p = reshape([p{:}], sz);
                return
            end

            state = char(state);
            if ~isscalar(state) && length(state) ~= length(qubits)
                error(message("quantum:QuantumMeasurement:invalidBitstring"));
            end

            if isscalar(state)
                state = repmat(state, 1, length(qubits));
            end

            states = obj.MeasuredStatesChar;
            probs = obj.Probabilities;

            % Reduce number of qubits to the requested qubits set
            states = states(:, qubits);
            [states, ~, id] = unique(states, 'rows');
            probs = accumarray(id, probs);

            ind = find(string(states) == state);

            if isempty(ind)
                if isempty(probs)
                    p = nan(1, 'like', probs);
                else
                    p = zeros(1, 'like', probs);
                end
            else
                p = probs(ind);
            end
        end
    end

    methods (Access=private)
        function thresh = defaultThreshold(obj)
            % Multiply by NumQubits assuming that the number of gates
            % applied to a state is roughly proportional to the number of
            % qubits.
            thresh = max(1, obj.NumQubits)*eps;
        end
    end

    properties(Constant, Access='private')
        % Version of the serialization and deserialization
        % format. This is used for managing forward compatibility. Value is
        % saved in 'versionSavedFrom' when an instance is serialized.
        %
        %   1.0 : original shipping version
        %   2.0 : add Probabilities property
        %   3.0 : support array
        version = 3.0;
    end
    methods (Hidden)
        function s = saveobj(m)
            % This is valid for the array case but only sees a scalar instance.
            if any(isnan(m.Counts))
                minCompatibleVersion = 2;
            else
                minCompatibleVersion = 1;
            end
            s = struct('MeasuredStates', m.MeasuredStates, ...
                'Counts', m.Counts, ...
                'Probabilities', m.Probabilities, ...
                'versionSavedFrom', quantum.gate.QuantumMeasurement.version, ...
                'minCompatibleVersion', minCompatibleVersion);
        end
    end
    methods(Hidden, Static)
        function m = loadobj(s)
            % This is valid for the array case but only sees a scalar instance.
            if quantum.gate.QuantumMeasurement.version < s.minCompatibleVersion
                id = 'quantum:QuantumMeasurement:IncompatibleVersion';
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.QuantumMeasurement', getString(message(id))));
                warning(id, loadWarningString);

                m = quantum.gate.QuantumMeasurement('', zeros(0, 1));
                return
            end

            try
                if s.versionSavedFrom == 1 || ~any(isnan(s.Counts))
                    m = quantum.gate.QuantumMeasurement(s.MeasuredStates, s.Counts);
                else
                    m = quantum.gate.QuantumMeasurement(s.MeasuredStates, s.Probabilities, "probs");
                end
            catch err
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'quantum.gate.QuantumMeasurement', err.message));
                warning(err.identifier, loadWarningString);

                m = quantum.gate.QuantumMeasurement('', zeros(0, 1));
                return
            end
        end
    end
end

function mustBeScalarQuantumMeasurement(obj)
if ~isa(obj, 'quantum.gate.QuantumMeasurement') || ~isscalar(obj)
    error(message("quantum:QuantumMeasurement:mustBeScalar"))
end
end

function mustBeStringVectorOrCharMatrix(states)
if ~(isstring(states) && isvector(states)) && ~(ischar(states) && ismatrix(states))
    error(message('quantum:QuantumMeasurement:invalidStates'))
end
end

function mustBeBasisString(str)
% Passing QuantumState instead of just numQubits as an FAV workaround
mustBeTextScalar(str)
str = char(str);
if ~all(ismember(str, '01'))
    error(message("quantum:QuantumMeasurement:invalidBitstring"));
end
end

function mustBeQubits(qubits, meas)
smallestQubit = min([meas.NumQubits]);
if ~isnumeric(qubits) || ~isvector(qubits) || ~isreal(qubits) || any(floor(qubits) ~= qubits, "all") || ...
        numel(unique(qubits)) ~= numel(qubits) || any(qubits < 1, "all")
    error(message("quantum:QuantumMeasurement:invalidQubits"));
end
if ~isempty(smallestQubit) && any(qubits > smallestQubit)
    error(message("quantum:QuantumMeasurement:invalidQubits"));
end
end
