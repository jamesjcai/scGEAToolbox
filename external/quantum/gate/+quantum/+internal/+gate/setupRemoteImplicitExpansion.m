function [bIdx, perm, szOutIdxPerm] = setupRemoteImplicitExpansion(szA, szB)
% Internal use only. Returns the expanded indices into two arrays of size 
% szA and szB. The bIdx is a matrix of size [numel(A) expandDimsB]. 
% The perm and szOutIdxPerm are permutations to undo the expansion.

% This is used for implicit expansion of run(c, dev, Observable=o) when
% c and o are non-scalar. It follows expansion rules to match observe(c,o).

%   Copyright 2025 The MathWorks, Inc.

[szA, szB, szOut] = quantum.internal.gate.setupImplicitExpansion(szA, szB);

% Index array for implicit expansion
bIdx = reshape(1:prod(szB), szB);

% Expand B to the same shape as the output
repInd = szOut;
repInd(szOut == szB) = 1;
bIdx = repmat(bIdx, repInd);

% Find a permutation that places all dimensions of A that are not 1 on the
% left, and all its scalar dimensions on the right. The szA, szB, and szOut
% all have same length.
scalarDimsA = szA==1;
perm = [find(~scalarDimsA) find(scalarDimsA)];

% Permute and reshape
bIdx = permute(bIdx, perm);

szOutIdxPerm = size(bIdx);

expandDimsB = prod(szB(scalarDimsA));
expandDimsA = prod(szA(~scalarDimsA));

bIdx = reshape(bIdx, expandDimsA, expandDimsB);

end