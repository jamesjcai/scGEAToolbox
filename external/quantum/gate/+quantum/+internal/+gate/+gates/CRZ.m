classdef CRZ < quantum.internal.gate.gates.ControlledQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        controlQubit
        targetQubit
        gateMatrix
        theta
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "RZ";
        type = "crz";
    end

    methods
        function obj = CRZ(controlQubit, targetQubit, theta)
            obj@quantum.internal.gate.gates.ControlledQubitGate(controlQubit, targetQubit);
            obj.gateMatrix = [exp(-1i*theta/2), 0;...
                0, exp(1i*theta/2)];
            obj.theta = theta;
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.CRZ(obj.controlQubit, obj.targetQubit, -obj.theta);
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
                def = [rzGate(2,obj.theta/2) cxGate(1,2) rzGate(2,-obj.theta/2) cxGate(1,2)];
                cg = compositeGate(def, [obj.controlQubit obj.targetQubit], Name=obj.type);
                [instr, def, map] = getQASM(cg, Unpack, map, SupportedGates, varargin{:});

                true_instr = quantum.internal.gate.formatQASMInstruction(obj);
                instr = "//start "+true_instr+...
                        instr+...
                        "//end "+true_instr;
            end
        end

        function a = angles(obj)
            a = obj.theta;
        end

        function obj = setAngles(obj, theta)
            obj.theta = theta;
            obj.gateMatrix = [exp(-1i*theta/2), 0;...
                0, exp(1i*theta/2)];
        end
    end

    methods(Hidden, Static)
        function [Nctrl, Ntrgt, Nangle] = numInputs
            Nctrl = 1;
            Ntrgt = 1;
            Nangle = 1;
        end
    end
end
