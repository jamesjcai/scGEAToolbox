function [aIdx, bIdx, szOut] = setupLocalImplicitExpansion(szA, szB)
% Internal use only. Returns the expanded indices into two arrays of size 
% szA and szB. The output aIdx and bIdx have the same size as szOut.
% 
% This is used for implicit expansion of simulate(c,s) and observe(c,o)
% when c,s,o are non-scalar. It follows expansion rules from the plus
% operator.
%
%   Copyright 2025 The MathWorks, Inc.

[szA, szB, szOut] = quantum.internal.gate.setupImplicitExpansion(szA, szB);

aIdx = reshape(1:prod(szA), szA);
bIdx = reshape(1:prod(szB), szB);
% Expand A indices
repIdx = szOut;
repIdx(szOut==szA) = 1;
aIdx = repmat(aIdx, repIdx);
% Expand B indices
repIdx = szOut;
repIdx(szOut==szB) = 1;
bIdx = repmat(bIdx, repIdx);

end