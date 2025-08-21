function [obj, maxStallIter] = initializeRunTimeDefaults(obj, numVar)
%INITIALIZERUNTIMEDEFAULTS Initialize run time default options
%
%   [OBJ, MAXSTALLITER] = INITIALIZERUNTIMEDEFAULTS(OBJ, NUMVAR)
%   initializes run time defaults for TabuSearch, given the number of
%   problem variables.

%   Copyright 2023 The MathWorks, Inc.

% Set maximum number of iterations
if isempty(obj.MaxIterations)
    obj.MaxIterations = obj.Mu*numVar;
end

if isempty(obj.MaxStallIterations)
    if numVar <= 1000
        maxStallIter = 0.01*obj.MaxIterations;
    else
        maxStallIter = obj.MaxIterations;
    end
else
    maxStallIter = obj.MaxStallIterations;
end



