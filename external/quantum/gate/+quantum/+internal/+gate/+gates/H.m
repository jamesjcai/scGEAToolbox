classdef H < quantum.internal.gate.gates.OneQubitGate
    %

    %   Copyright 2021-2022 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix = 1/sqrt(2)*[1 1; 1 -1];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "H";
        type = "h";
    end

    methods
        function obj = H(targetQubit)
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
