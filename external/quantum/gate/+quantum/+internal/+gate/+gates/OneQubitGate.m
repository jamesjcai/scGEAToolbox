classdef OneQubitGate < quantum.internal.gate.gates.Gate

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(Abstract, SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix
    end

    properties(Abstract, SetAccess = private, Hidden, GetAccess = public)
        type
    end

    methods
        % constructor
        function obj = OneQubitGate(targetQubit)
            obj.targetQubit = targetQubit;
        end
    end

    methods(Hidden)
        function outputState = applyToState(obj, inputState, numQubits)
            outputState = quantum.internal.gate.applyMatToDim(inputState,...
                obj.gateMatrix, numQubits+1-obj.targetQubit);
        end

        function [instr, def, map] = getQASM(obj, ~, map, varargin)
            instr = quantum.internal.gate.formatQASMInstruction(obj);
            def = ""; 
        end

        function qb = getQubits(obj)
            qb = obj.targetQubit;
        end

        function obj = setQubits(obj, qb)
            assert(isscalar(qb));
            obj.targetQubit = qb;
        end

        function qb = getTargetQubits(obj)
            qb = obj.targetQubit;
        end

        function qb = getControlQubits(obj)
            qb = zeros(1, 0);
        end
    end
end
