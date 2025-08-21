function invobj = inv(obj)
%INV  Inverse of a quantumCircuit.
%
%   invcirc = INV(circ) returns a new quantumCircuit invcirc with inverse
%   behavior of the original quantumCircuit circ. The order of the returned
%   gates is reversed from the original gates and each original gate is
%   replaced by its inverse.
%
%   Example:
%       % Construct and invert quantumCircuit.
%       gates = [hGate(1); sGate(2); rxGate(1, pi/3)];
%       circ = quantumCircuit(gates);
%       invCirc = inv(circ);
%       circ.Gates
%       invCirc.Gates
%
%   See also quantumCircuit, quantumCircuit/getMatrix,
%   quantum.gate.SimpleGate/inv, quantum.gate.CompositeGate/inv

%   Copyright 2021-2022 The MathWorks, Inc.

cG = obj.Gates;
cG = flip(cG);
for ii=1:length(cG)
    cG(ii) = inv(cG(ii));
end
invobj = quantumCircuit(obj.NumQubits);
invobj.Gates = cG;

% Adjust the name of the quantumCircuit by adding / removing "_inv".
if strlength(obj.Name) ~= 0
    if endsWith(obj.Name, "_inv")
        n = char(obj.Name);
        n(end-3:end) = [];
        invobj.Name = string(n);
    else
        invobj.Name = obj.Name + "_inv";
    end
end

if nargout == 0
    warning(message('quantum:quantumCircuit:NotModified'));
end
end