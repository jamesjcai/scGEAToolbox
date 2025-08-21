function checkQubitBlocksForPlot(numQubits, args)
%

% Used in quantumCircuit/plot and CompositeGate/plot for checking the
% QubitBlocks NVP

%   Copyright 2023 The MathWorks, Inc.

if rem(numel(args), 2) ~= 0
    % This will error in the constructor
    return
end

for ii=1:2:numel(args)
    name = args{ii};
    value = args{ii+1};
    % Check at least 6 characters to avoid matching ambiguous input 'Qubit'
    % which could be either 'QubitBlocks' or 'QubitLabelLocation'.
    if strncmpi(name, 'QubitBlocks', max(strlength(name), 6))
        quantum.internal.gate.checkQubitBlocks(value, numQubits);
    end
end
end