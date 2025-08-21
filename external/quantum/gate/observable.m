classdef observable
    %OBSERVABLE Observable of a quantumCircuit
    %
    %   An observable represents a set of projective qubit measurements,
    %   expressed as weighted Pauli strings or a Hermitian matrix. The
    %   Pauli strings specify the measurement basis on each qubit. They are
    %   equivalent to the Hermitian representation, because each Pauli
    %   string represents the Kronecker product of the corresponding Pauli
    %   matrices, such that the weighted sum is the Hermitian matrix.
    %
    %   obs = observable(paulis, weights) takes a string array or char matrix
    %   paulis and numeric array weights of the same length. All elements in
    %   the paulis input must be "X", "Y", "Z", or "I" and all numbers in
    %   input weights must be real.
    %
    %   obs = observable(matrix) takes a Hermitian matrix of size
    %   2^N by 2^N, where N is the number of qubits. The Pauli strings and
    %   weights are computed from the input matrix.
    %
    %   Example:
    %       % Construct an observable for 2 qubits, where X is observed on
    %       % the 1st and 2nd qubit with weight 2.4, and Z is observed on
    %       % the first qubit with weight -1.2. Get the Hermitian matrix of
    %       % this observable.
    %       obs = observable(["XX" "ZI"], [2.4 -1.2]);
    %       H = getMatrix(obs);
    %
    %   Example:
    %       % Construct an observable for 1 qubit described by a Hermitian
    %       % matrix.
    %       H = [2.5 0.3-1.2j; 0.3+1.2j -2.3]
    %       obs = observable(H);
    %
    %   See also observe, qubo2ising

    %   Copyright 2025 The MathWorks, Inc.
    properties
        Paulis
        Weights
        NumQubits
    end

    properties(Access=private)
        Matrix
    end

    properties(Hidden)
        IsDiagonal
    end

    methods

        function obj = observable(varargin)
            narginchk(1,2)
            if nargin==1
                % Matrix Constructor
                H = varargin{1};
                [paulis, weights, N, isDiagonal] = matrix2pauli(H);
                obj.Matrix = H;
            else
                % Pauli Constructor
                paulis = varargin{1};
                weights = varargin{2};
                [paulis, weights, N, isDiagonal] = checkPauli(paulis, weights);
            end

            obj.Paulis = paulis;
            obj.Weights = weights;
            obj.NumQubits = N;
            obj.IsDiagonal = isDiagonal;
        end

        function H = getMatrix(obj)
            % Returns sparse matrix

            if ~isempty(obj.Matrix)
                H = sparse(obj.Matrix);
                return
            end

            % Empty cases
            N = obj.NumQubits;
            if isempty(obj.Paulis)
                % Only possible when Weights are 0
                H = sparse(2^N, 2^N);
                return
            elseif isequal(obj.Paulis, "")
                % Only possible when N is 0
                H = sparse(obj.Weights);
                return
            end

            % Index of each Pauli character into map
            if obj.IsDiagonal
                map = 'IZ';
            else
                map = 'IXYZ';
            end
            [~, mapIdx] = ismember(char(obj.Paulis), map);

            % Index of each term into the reshaped matrix
            if N==1
                linIdx = mapIdx;
            else
                mapIdx = num2cell(fliplr(mapIdx), 1);
                linIdx = sub2ind(length(map)*ones(1,N), mapIdx{:});
            end

            if obj.IsDiagonal
                % Project in full for performance
                twoU = [1 1; 1 -1];
                weights = zeros(2^N, 1, 'like', obj.Weights);
                weights(linIdx) = obj.Weights;
                h = fullProjDims(weights, twoU, N);
                H = spdiags(h, 0, 2^N, 2^N);
            else
                twoU = sparse([1 0 0 1; 0 1 1j 0; 0 1 -1j 0; 1 0 0 -1]);
                [linIdx, h] = sparseProjDims(linIdx, obj.Weights, twoU, N);
                H = sparseMortonInvPermute(linIdx, h, N);
            end
        end
    end

    properties(Constant, Access=private)
        % Version of the serialization and deserialization
        % format. This is used for managing forward compatibility. Value is
        % saved in 'versionSavedFrom' when an instance is serialized.
        %
        %   1.0 : original shipping version
        version = 1.0;
    end
    methods(Hidden)
        function s = saveobj(obs)
            % Workaround to ensure proper loading for empty cases
            if isequal(obs.Paulis, "")
                paulis = obs.Paulis;
            else
                paulis = reshape(char(obs.Paulis), [], obs.NumQubits);
            end
            s = struct('Paulis', paulis, ...
                'Weights', obs.Weights, ...
                'versionSavedFrom', observable.version, ...
                'minCompatibleVersion', 1);
        end

        function tf = isequal(obs, other)
            % Overload so both constructors return equal objects when their
            % visible properties match. 
            if ~isa(other, 'observable')
                tf = false;
            else
                tf = isequal(obs.NumQubits, other.NumQubits) && ...
                     isequal(obs.Paulis, other.Paulis) && ...
                     isequal(obs.Weights, other.Weights);
            end
        end
    end
    methods(Hidden, Static)
        function obs = loadobj(s)
            if observable.version < s.minCompatibleVersion
                id = 'quantum:observable:IncompatibleVersion';
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'observable', getString(message(id))));
                warning(id, loadWarningString);

                obs = observable(0);
                return
            end

            try
                obs = observable(s.Paulis, s.Weights);
            catch err
                loadWarningString = getString(message('MATLAB:load:classError', ...
                    'observable', err.message));
                warning(err.identifier, loadWarningString);

                obs = observable(0);
                return
            end
        end
    end
end

%% Helper Functions

function [pauli, weights, N, isDiagonal] = matrix2pauli(H)
% Project input matrix to coefficients in the Pauli basis

if isempty(H) || ~ishermitian(H)
    error(message("quantum:observable:invalidHermitian"))
end

if isscalar(H)
    N = 0;
    [~, weights] = removeSmallWeights(H, N);
    if isempty(weights)
        pauli = string.empty(0,1);
    else
        pauli = "";
    end
    isDiagonal = true;
    return
end

N = log2(size(H,1));
if floor(N)~=ceil(N)
    error(message("quantum:observable:invalidHermitian"))
end

isDiagonal = isdiag(H);

% Reshape and project based on type for performance
if isDiagonal
    U = [0.5 0.5; 0.5 -0.5];
    map = 'IZ';
    h = full(diag(H));
    weights = fullProjDims(h, U, N);
    [linIdx, weights] = removeSmallWeights(weights, N);
else
    % U = reshape(cat(3, I, X, Y, Z), 4, 4)/2
    U = [0.5 0 0 0.5; 0 0.5 0.5 0; 0 -0.5j 0.5j 0; 0.5 0 0 -0.5];
    map = 'IXYZ';
    if issparse(H)
        U = sparse(conj(U));
        [linIdx, h] = sparseMortonPermute(H, N);
        [hIdx, weights] = sparseProjDims(linIdx, h, U, N);
        [wIdx, weights] = removeSmallWeights(weights, N);
        linIdx = hIdx(wIdx);
    else
        h = fullMortonPermute(H, N);
        weights = fullProjDims(h, U, N);
        [linIdx, weights] = removeSmallWeights(weights, N);
    end
end

% Build Pauli strings
if N==1
    mapIdx = linIdx;
    pauli = string(map(mapIdx).');
else
    cc = cell(1, N);
    [cc{:}] = ind2sub(length(map)*ones(1, N), linIdx);
    mapIdx = horzcat(cc{end:-1:1});
    pauli = string(map(mapIdx));
end
end


function [p, w, N, isDiagonal] = checkPauli(p, w)

% Check datatypes
if ~(ischar(p) || isstring(p))
    error(message("quantum:observable:invalidPauli"))
end
if ~(isfloat(w) && isreal(w) && allfinite(w) && isvector(w))
    error(message("quantum:observable:invalidWeights"))
end

if isstring(p)
    if ~isvector(p) || any(ismissing(p))
        error(message("quantum:observable:invalidPauliString"))
    end

    numTerms = length(p);

    % Length of each string element must be the same
    N = strlength(p);
    if isempty(N)
        N = 0;
    else
        if ~all(N == N(1))
            error(message("quantum:observable:invalidPauliString"))
        end
        N = N(1);
    end

    % Convert to char for input checking below
    p = char(reshape(p, [], 1));
else
    if ~ismatrix(p)
        error(message("quantum:observable:invalidPauliChar"))
    end
    numTerms = size(p,1);
    N = size(p, 2);
end

w = reshape(w, [], 1);
if numTerms~=length(w)
    error(message("quantum:observable:invalidWeights"))
end

if N==0
    isDiagonal = true;
    p = strings(numTerms, 1);
else
    % Verify each character of the char matrix
    p = upper(p);
    [tf, idx] = ismember(p, 'IXYZ');
    if ~all(tf, "all")
        error(message("quantum:observable:invalidPauli"))
    end
    % Combinations of "I" and "Z" are always diagonal
    isDiagonal = all(ismember(idx, [1 4]), "all");
    p = string(p);
end

% Pauli strings are unique and sorted alphabetically
if ~issorted(p, 'strictascend')
    [p, ~, idx] = unique(p);
    w = accumarray(idx, w);
end

[linIdx, w] = removeSmallWeights(w, N);
p = p(linIdx);
end

function h = fullMortonPermute(H, N)
% Reshape input matrix to vector using the Morton (Z-order) pattern. 
h = ipermute(reshape(H, 2*ones(1, 2*N)), [1:2:2*N 2:2:2*N]);
end

function [linIdx, h] = sparseMortonPermute(H, N)
% Reshape input matrix to vector using the Morton (Z-order) pattern. 
% Sparse implementation determines linear indices by interleaving bits of
% the rows and columns.
[ii, jj, h] = find(H);
linBin= repmat('0', [length(h) 2*N]);
linBin(:, 1:2:end) = dec2bin(ii-1, N);
linBin(:, 2:2:end) = dec2bin(jj-1, N);
linIdx = bin2dec(linBin)+1;
end

function H = sparseMortonInvPermute(linIdx, h, N)
% Reshape input vector to matrix using the inverse Morton (Z-order) pattern.  
llbin = dec2bin(linIdx-1, 2*N);
rows = bin2dec(llbin(:, 2:2:end))+1;
cols = bin2dec(llbin(:, 1:2:end))+1;
H = sparse(rows, cols, h, 2^N, 2^N);
end

function h = fullProjDims(h, U, N)

sz = size(U,1);
h = U*reshape(h, sz, []); % dim=1
for dim=2:N
    h = reshape(h, sz^(dim-1), sz, []);
    h = pagemtimes(h, U.');
end
h = h(:);
end

function [ind, h] = sparseProjDims(ind, h, U, N)

sz = 4;
p = ceil(N/2);
q = N - p;

% Reshape H to be a close-to-square matrix
H = sparse(ind, 1, h, sz^N, 1);
H = reshape(H, sz^p, sz^q);

% Apply U to all dimensions in the rows and columns of H
id = @(n) speye(sz^n);
for ii=1:p
    Ukron = kron(kron(id(ii-1), U), id(p-ii));
    H = Ukron * H;
end
for ii=1:q
    Ukron = kron(kron(id(ii-1), U), id(q-ii));
    H = H * Ukron.';
end

H = reshape(H, [], 1);
[ind, ~, h] = find(H);
end

function [idx, w] = removeSmallWeights(w, N)
% Remove weights close to 0.
thres = defaultTolerance(w, N);
idx = find(abs(w)>thres);
% Ensure [] becomes empty(0,1)
idx = idx(:);
w = w(idx);
w = w(:);
end

function tol = defaultTolerance(values, numQubits)
% Tolerance for removing small weights in constructors
tol = max(abs(values))*sqrt(2^numQubits*eps(class(values)));
end
