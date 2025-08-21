classdef X < quantum.internal.gate.gates.OneQubitGate
    %

    %   Copyright 2021-2022 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix = [0 1; 1 0];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "X";
        type = "x";
    end

    methods
        function obj = X(targetQubit)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
        end

        function obj = inv(obj)
        end
    end

    methods(Hidden, Static)
        function [Nctrl, Ntrgt, Nangle] = numInputs
            Nctrl = 0;
            Ntrgt = 1;
            Nangle = 0;
        end
    end
end
