function out = undoRemoteImplicitExpansion(out, perm, szOutIdxPerm)
% Internal use only. Returns the array out with the expected size from
% implicit expansion. This undoes the permutation done by setupRemoteImplicitExpansion
%
%   Copyright 2025 The MathWorks, Inc.
out = reshape(out, szOutIdxPerm);
out = ipermute(out, perm);
end