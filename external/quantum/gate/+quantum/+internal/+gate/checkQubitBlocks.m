function checkQubitBlocks(qubitBlocks, numQubits)
%

% Used in for checking the QubitBlocks property / NVP. Second input only
% used for the NVP to plot; when setting the property of
% QuantumCircuitChart to have invalid sum, only a warning is given.

%   Copyright 2023 The MathWorks, Inc.

if ~isvector(qubitBlocks) || ~isnumeric(qubitBlocks) || ...
         ~isreal(qubitBlocks) || ~allfinite(qubitBlocks) || ...
        ~all(floor(qubitBlocks) == qubitBlocks) || ~all(qubitBlocks >= 0)
    error(message('quantum:quantumCircuit:plotInvalidQubitBlocks'))
end

if nargin > 1 && sum(qubitBlocks) ~= numQubits
    error(message('quantum:quantumCircuit:SumQubitBlocks'))
end
