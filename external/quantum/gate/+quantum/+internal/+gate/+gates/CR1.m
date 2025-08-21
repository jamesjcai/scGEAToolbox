classdef CR1 < quantum.internal.gate.gates.ControlledQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        controlQubit
        targetQubit
        gateMatrix
        lambda
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "R1";
        type = "cr1";
    end

    methods
        function obj = CR1(controlQubit, targetQubit, lambda)
            obj@quantum.internal.gate.gates.ControlledQubitGate(controlQubit, targetQubit);
            obj.gateMatrix = [1 0; 0 exp(1i*lambda)];
            obj.lambda = lambda;
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.CR1(obj.controlQubit, obj.targetQubit, -obj.lambda);
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
                % NOTE: Adds global phase exp(1j*-lambda/4)
                gates = [crzGate(1, 2, obj.lambda) rzGate(1, obj.lambda/2)];
                cg = compositeGate(gates, [obj.controlQubit obj.targetQubit], Name=obj.type);
                [instr, def, map] = getQASM(cg, Unpack, map, SupportedGates, varargin{:});

                true_instr = quantum.internal.gate.formatQASMInstruction(obj);
                instr = "//start "+true_instr+...
                        instr+...
                        "//end "+true_instr;
            end
        end
        
        function a = angles(obj)
            a = obj.lambda;
        end

        function obj = setAngles(obj, lambda)
            obj.lambda = lambda;
            obj.gateMatrix = [1 0; 0 exp(1i*lambda)];
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
