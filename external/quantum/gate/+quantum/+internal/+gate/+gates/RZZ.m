classdef RZZ < quantum.internal.gate.gates.TwoQubitRotationGate
    %

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit1
        targetQubit2
        gateMatrix
        theta
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "RZZ";
        type = "rzz";
    end

    methods
        function obj = RZZ(targetQubit1, targetQubit2, theta)
            obj@quantum.internal.gate.gates.TwoQubitRotationGate(targetQubit1, targetQubit2);
            obj.theta = theta;
            obj.gateMatrix = diag([exp(-1i*theta/2) exp(1i*theta/2) exp(1i*theta/2) exp(-1i*theta/2)]);
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.RZZ(obj.targetQubit1, obj.targetQubit2, -obj.theta);
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
                gates = [cxGate(1,2) rzGate(2, obj.theta) cxGate(1,2)];
                cg = compositeGate(gates, [obj.targetQubit1 obj.targetQubit2], Name=obj.type);
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
            obj.gateMatrix = diag([exp(-1i*theta/2) exp(1i*theta/2) exp(1i*theta/2) exp(-1i*theta/2)]);
        end
    end

    methods(Hidden, Static)
        function [Nctrl, Ntrgt, Nangle] = numInputs
            Nctrl = 0;
            Ntrgt = 2;
            Nangle = 1;
        end
    end
end