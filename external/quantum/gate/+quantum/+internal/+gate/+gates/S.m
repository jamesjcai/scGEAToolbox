classdef S < quantum.internal.gate.gates.OneQubitGate
    %

    %   Copyright 2021-2022 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix = [1 0;0 1i];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "S";
        type = "s";
    end

    methods
        function obj = S(targetQubit)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.SI(obj.targetQubit);
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

