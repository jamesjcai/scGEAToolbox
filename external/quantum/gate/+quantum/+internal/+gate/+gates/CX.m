classdef CX < quantum.internal.gate.gates.ControlledQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        controlQubit
        targetQubit
        gateMatrix = [0 1; 1 0];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "X";
        type = "cx";
    end

    methods
        function obj = CX(controlQubit, targetQubit)
            obj@quantum.internal.gate.gates.ControlledQubitGate(controlQubit, targetQubit);
        end

        function obj = inv(obj)
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
            Nctrl = 1;
            Ntrgt = 1;
            Nangle = 0;
        end
    end
end
