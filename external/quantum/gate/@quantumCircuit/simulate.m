function outputState = simulate(obj,inputState)
%SIMULATE  Simulate a quantumCircuit
%
%   S = SIMULATE(circ) simulates the quantum circuit circ and returns the
%   final state S after running the circuit as a QuantumState object. Each
%   qubit is initially in the |0> state.
%
%   S = SIMULATE(circ, inputState) additionally specifies the initial
%   state inputState passed to the circuit circ. Input state can be a
%   QuantumState object, basis string, or numeric vector. The default
%   input state is "0...0" with as many "0" characters as qubits in the
%   circuit.
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

%   Copyright 2021-2024 The MathWorks, Inc.

if isinf(2^obj.NumQubits)
    % Give useful error in case of overflow
    error(message('MATLAB:pmaxsize'));
end

if nargin < 2
    % Defaults to all |0> if no input state is specified
    inputState = zeros(2^obj.NumQubits,1);
    inputState(1) = 1;
else
    if ~isa(inputState, 'quantum.gate.QuantumState')
        inputState = quantum.gate.QuantumState(inputState);
    end

    if inputState.NumQubits~=obj.NumQubits
        error(message("quantum:quantumCircuit:simulateIncorrectNumQubits"))
    end
    inputState = inputState.Amplitudes;
end

outputState = inputState;
for g = 1:length(obj.Gates)
    outputState = applyToState(obj.Gates(g), outputState, obj.NumQubits);
end

outputState = quantum.gate.QuantumState(outputState(:));
end
