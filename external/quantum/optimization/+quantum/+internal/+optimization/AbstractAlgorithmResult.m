classdef AbstractAlgorithmResult
%ABSTRACTALGORITHMRESULT Base class for AlgorithmResult classes

%   Copyright 2022-2023 The MathWorks, Inc.    
    properties (SetAccess = protected, GetAccess = public)
        ClockTime
    end
    methods(Abstract=true)
        obj = update(obj)
        x = retrieveBestPoint(obj)
        obj = adjustForConstant(obj, constant)
    end
end