function mustBeDouble(input)
%MUSTBEDOUBLE Validate that value is double
%
%   MUSTBEDOUBLE(INPUT) throws an error if INPUT contains non-double
%   values. 

%   Copyright 2025 The MathWorks, Inc.

if ~isa(input, "double")
    error(message("quantum:annealing:validators:mustBeDouble"));
end
end
