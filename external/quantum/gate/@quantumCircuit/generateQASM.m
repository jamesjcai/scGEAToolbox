function qasm = generateQASM(obj, NameValueArgs)
%GENERATEQASM  Generate OpenQASM code for a quantumCircuit.
%
%   qasm = generateQASM(circ) returns a string representing OpenQASM code
%   for quantumCircuit circ.
%
%   Example:
%       % Construct and generate OpenQASM code for quantumCircuit.
%       gates = [hGate(1); cxGate(1, 2)];
%       circ = quantumCircuit(gates);
%       str = generateQASM(circ)
%
%   See also quantumCircuit, quantumCircuit/run

%   Copyright 2021-2025 The MathWorks, Inc.

arguments
    obj quantumCircuit
    NameValueArgs.Unpack {mustBeA(NameValueArgs.Unpack, 'logical')} = false
    NameValueArgs.SupportedGates {mustBeA(NameValueArgs.SupportedGates, 'dictionary')} 
end

if ~isscalar(obj)
    error(message('quantum:quantumCircuit:mustBeScalar'))
end

if ~isfield(NameValueArgs, 'SupportedGates')

    % Supported subset of the OpenQASM 3.0 stdgates.inc library
    SupportedGates = dictionary(["id","x","y","z","h","s","si","t","ti","r1","rx","ry","rz","cx","cy","cz","crx","cry","crz","ch","swap","ccx"], ...
                                ["id","x","y","z","h","s","sdg","t","tdg","u1","rx","ry","rz","cx","cy","cz","crx","cry","crz","ch","swap","ccx"]);
else

% Dictionary mapping MATLAB names to vendor names
SupportedGates = NameValueArgs.SupportedGates;

% Gates not supported by a vendor are expressed in terms of gates from
% the following set. If any of these gates are also unsupported, valid
% code cannot be generated.
minGateSet=["cx","h","rx","ry","rz","s","si","swap","t","ti","x","y","z"];

isMinGateSupported = isKey(SupportedGates, minGateSet);
if ~all(isMinGateSupported)
    error(message('quantum:generateQASM:missingSupportedGates'))
end
end

instructions = "";
definitions = "";
DefinedGates = containers.Map;

for i = 1:length(obj.Gates)
    [instr, def, DefinedGates] = getQASM(obj.Gates(i), NameValueArgs.Unpack, DefinedGates, SupportedGates); 
    instructions = instructions + instr;
    definitions = definitions + def;
end

header = "OPENQASM 3.0;"+newline+...
         'include "stdgates.inc";'+newline;

qreg = newline+sprintf("qubit[%g] q;", obj.NumQubits);
creg = newline+sprintf("bit[%g] c;", obj.NumQubits)+newline+newline;

measurements = "c = measure q;";

qasm = header+...
    definitions+...
    qreg+creg+...
    instructions+...
    measurements;
end




