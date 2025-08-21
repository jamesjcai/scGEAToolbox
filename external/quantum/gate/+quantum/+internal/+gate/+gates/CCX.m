classdef CCX < quantum.internal.gate.gates.ControlledControlledQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        controlQubit1
        controlQubit2
        targetQubit
        gateMatrix = [0 1; 1 0];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "X";
        type = "ccx";
    end

    methods
        function obj = CCX(controlQubit1, controlQubit2, targetQubit)
            obj@quantum.internal.gate.gates.ControlledControlledQubitGate(...
                controlQubit1, controlQubit2, targetQubit);
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
                gates = [hGate(3)
                    cxGate(2, 3)
                    tiGate(3)
                    cxGate(1, 3)
                    tGate(3)
                    cxGate(2, 3)
                    tiGate(3)
                    cxGate(1, 3)
                    tGate(2)
                    tGate(3)
                    cxGate(1, 2)
                    hGate(3)
                    tGate(1)
                    tiGate(2)
                    cxGate(1, 2)];
                cg = compositeGate(gates, [obj.controlQubit1 obj.controlQubit2 obj.targetQubit], Name=obj.type);
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
            Nctrl = 2;
            Ntrgt = 1;
            Nangle = 0;
        end
    end
end