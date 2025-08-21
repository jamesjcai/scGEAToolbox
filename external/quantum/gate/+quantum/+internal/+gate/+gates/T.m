classdef T < quantum.internal.gate.gates.OneQubitGate
    %

    %   Copyright 2021-2022 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix = [1 0;0 exp(1i*pi/4)];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "T";
        type = "t";
    end

    methods
        function obj = T(targetQubit)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.TI(obj.targetQubit);
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