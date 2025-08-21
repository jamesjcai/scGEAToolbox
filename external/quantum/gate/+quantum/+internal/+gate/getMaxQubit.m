function maxQubit = getMaxQubit(gates)
%

%   Copyright 2021-2022 The MathWorks, Inc.

maxQubit = 0;
for ii=1:length(gates)
    qb = getQubits(gates(ii));
    if ~isempty(qb)
        maxQubit = max(maxQubit, max(qb));
    end
end
end
