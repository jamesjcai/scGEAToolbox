classdef ID < quantum.internal.gate.gates.OneQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix = [1 0; 0 1];
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "ID";
        type = "id";
    end

    methods
        function obj = ID(targetQubit)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
        end

        function obj = inv(obj)
        end
    end

    methods(Hidden)
        function [instr, def, map] = getQASM(obj, ~, map, SupportedGates, varargin)
            
            def = "";

            if isKey(SupportedGates, obj.type) 
                % Supported: Use name recognized by vendor
                instr = quantum.internal.gate.formatQASMInstruction(obj, SupportedGates(obj.type));
            else
                % Not Supported: 
                true_instr = quantum.internal.gate.formatQASMInstruction(obj);
                instr = "//start "+true_instr+...
                        "//do nothing"+newline+...
                        "//end "+true_instr;
            end
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
