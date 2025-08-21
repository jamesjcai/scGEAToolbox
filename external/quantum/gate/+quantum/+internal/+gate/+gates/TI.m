classdef TI < quantum.internal.gate.gates.OneQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix = [1 0;0 exp(-1i*pi/4)];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "TI";
        type = "ti";
    end

    methods
        function obj = TI(targetQubit)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.T(obj.targetQubit);
        end
    end
    
    methods(Hidden)
        function [instr, def, map] = getQASM(obj, ~, map, SupportedGates, varargin)
            
            % Assumed Supported: Use name recognized by vendor
            instr = quantum.internal.gate.formatQASMInstruction(obj, SupportedGates(obj.type));
            def = "";
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
