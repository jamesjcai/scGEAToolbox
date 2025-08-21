classdef CH < quantum.internal.gate.gates.ControlledQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.

    properties(SetAccess = protected, GetAccess = public)
        controlQubit
        targetQubit
        gateMatrix = 1/sqrt(2)*[1 1; 1 -1];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "H";
        type = "ch";
    end

    methods
        function obj = CH(controlQubit, targetQubit)
            obj@quantum.internal.gate.gates.ControlledQubitGate(controlQubit, targetQubit);
        end

        function obj = inv(obj)
        end
    end

    methods(Hidden)
        function [instr, def, map] = getQASM(obj, Unpack, map, SupportedGates, varargin)
                        
            if isKey(SupportedGates, obj.type) 
                % Supported: Use name recognized by vendor
                instr = quantum.internal.gate.formatQASMInstruction(obj, SupportedGates(obj.type));
                def = ""; 
            else
                % Not Supported: Implement as other gates 
                gates = [ryGate(2, pi/4) cxGate(1,2) ryGate(2, -pi/4)];
                cg = compositeGate(gates, [obj.controlQubit obj.targetQubit], Name=obj.type);
                [instr, def, map] = getQASM(cg, Unpack, map, SupportedGates, varargin{:});

                true_instr = quantum.internal.gate.formatQASMInstruction(obj);
                instr = "//start "+true_instr+...
                        instr+...
                        "//end "+true_instr;
            end
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
