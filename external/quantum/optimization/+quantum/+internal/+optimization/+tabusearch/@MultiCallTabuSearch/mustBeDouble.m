function mustBeDouble(input)
%MUSTBEDOUBLE Validate that value is double
%
%   MUSTBEDOUBLE(INPUT) throws an error if A is not a double.

%   Copyright 2023 The MathWorks, Inc.

if ~isa(input, "double")
    throwAsCaller(MException(message('quantum:annealing:TabuSearch:mustBeDouble')));
end

end
