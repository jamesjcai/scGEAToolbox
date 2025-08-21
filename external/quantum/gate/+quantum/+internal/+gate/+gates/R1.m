classdef R1 < quantum.internal.gate.gates.OneQubitGate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix
        lambda
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "R1";
        type = "r1";
    end

    methods
        function obj = R1(targetQubit,lambda)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
            obj.lambda = lambda;
            obj.gateMatrix = [1 0; 0 exp(1i*lambda)];
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.R1(obj.targetQubit, -obj.lambda);
        end
    end

    methods (Hidden)
        function [instr, def, map] = getQASM(obj, Unpack, map, SupportedGates, varargin)
                        
            if isKey(SupportedGates, obj.type) 
                % Supported: Use name recognized by vendor
                instr = quantum.internal.gate.formatQASMInstruction(obj, SupportedGates(obj.type));
                def = ""; 
            else
                % Not Supported: Implement as other gates
                % NOTE: Adds global phase exp(-1i*lambda) 
                cg = compositeGate(rzGate(1,obj.lambda), [obj.targetQubit], Name=obj.type);
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
            Nctrl = 0;
            Ntrgt = 1;
            Nangle = 1;
        end
    end
end
