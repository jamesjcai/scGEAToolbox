function mustHaveSameNumElements(v, w)
%MUSTHAVESAMENUMELEMENTS Validate that values have same number of elements
%
%   MUSTHAVESAMENUMELEMENTS(V, W) throws an error if V and W have different
%   number of elements.

%   Copyright 2025 The MathWorks, Inc.

numVars = numel(v);
if numel(w) ~= numVars
    error(message("quantum:annealing:validators:mustHaveSameNumElements"));        
end

end
