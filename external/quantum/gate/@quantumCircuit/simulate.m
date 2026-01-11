function outputState = simulate(obj,inputState)
%SIMULATE  Simulate a quantumCircuit
%
%   S = SIMULATE(circ) simulates the quantum circuit circ and returns the
%   final state S as a QuantumState object. Each qubit is initially in the
%   |0> state. The output state array S has the same size as the input
%   circuit array circ.
%
%   S = SIMULATE(circ, inputState) additionally specifies the initial
%   state inputState passed to the circuit circ. Input state can be a
%   QuantumState object, basis string, or numeric vector. The default
%   input state is "0...0" with as many "0" characters as qubits in the
%   circuit. The basis string and numeric vector input must represent the
%   same number of qubit for all circuits. The QuantumState input array
%   must match the number of qubits in each corresponding circuit.
%
%   Example:
%       % Construct and simulate quantumCircuit.
%       gates = [hGate(1); cxGate(1, 2)];
%       circ = quantumCircuit(gates);
%       s = simulate(circ)
%
%   Example:
%       % Construct quantumCircuit and simulate with various input states.
%       c = quantumCircuit(cxGate(1, 2));
%       for ket=["00", "01", "10", "11"]
%           s = simulate(c, ket);
%           disp("|" + ket + "> -> " + formula(s));
%       end
%
%   See also quantumCircuit, quantum.gate.QuantumState, quantumCircuit/run

%   Copyright 2021-2025 The MathWorks, Inc.

if isempty(obj)
    outputState = quantum.gate.QuantumState.empty(size(obj));
    return
end

if ~allfinite(2.^[obj.NumQubits])
    % Give useful error in case of overflow
    error(message('MATLAB:pmaxsize'));
end

% Determine inputState
if nargin < 2
    % Defaults to all |0> for each circuit when no input state is specified.
    L = numel(obj);
    szC = size(obj);
    inputState = cell(szC);
    N = [obj.NumQubits];
    if all(N(:)==N(1))
        % Only need 1 object
        inputState = quantum.gate.QuantumState(repmat('0',[1 N(1)]));
        sIdx = ones(1, L);
    else
        for ii = 1:L
            inputState{ii} = quantum.gate.QuantumState(repmat('0',[1 obj(ii).NumQubits]));
        end
        inputState = reshape([inputState{:}], szC);
        sIdx = 1:L;
    end
    cIdx = 1:L;
    % No expansion
    szOut = szC;
else
    if ~isa(inputState, 'quantum.gate.QuantumState')
        inputState = quantum.gate.QuantumState(inputState);
    end
    
    % Determine output size
    [cIdx, sIdx, szOut] = quantum.internal.gate.setupLocalImplicitExpansion(size(obj), size(inputState));
    
    L = prod(szOut);
    nC = [obj.NumQubits];
    nS = [inputState.NumQubits];
    n1 = nC(cIdx);
    n2 = nS(sIdx);
    if ~isequal(n1(:), n2(:))
        error(message("quantum:quantumCircuit:simulateIncorrectNumQubits"))
    end
end
   
if isscalar(obj)
    outputState = scalarsimulate(obj, inputState);
else
    outputState = cell(szOut);
    for ii = 1:L
        c = obj(cIdx(ii));
        s = inputState(sIdx(ii));
        outputState{ii} = scalarsimulate(c, s);
    end
    outputState = reshape([outputState{:}], szOut);
end

end


function xo = scalarsimulate(c, xi)
    g = c.Gates;
    xo = xi.Amplitudes;
    n = c.NumQubits;
    for jj = 1:length(g)
        xo = applyToState(g(jj), xo, n);
    end
    xo = quantum.gate.QuantumState(xo(:));
end