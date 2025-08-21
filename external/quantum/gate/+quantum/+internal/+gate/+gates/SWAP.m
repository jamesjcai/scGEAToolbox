classdef SWAP < quantum.internal.gate.gates.Gate
    %

    %   Copyright 2021-2022 The MathWorks, Inc.

    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        targetQubit2
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "SWAP";
        type = "swap";
    end

    methods
        %constructor
        function obj = SWAP(targetQubit, targetQubit2)
            if targetQubit == targetQubit2
                error(message('quantum:gates:matchingQubits'));
            end
            obj.targetQubit = targetQubit;
            obj.targetQubit2 = targetQubit2;
        end

        function obj = inv(obj)
        end
    end

    methods(Hidden)
        function outputState = applyToState(obj, inputState, numQubits)
            if numel(inputState) < 4
                error(message('quantum:gates:tooFewQubits', 2));
            end
%             cnot = CX(obj.targetQubit, obj.targetQubit2);
%             rev_cnot = CX(obj.targetQubit2, obj.targetQubit);
% 
%             outputState = applyToState(cnot, inputState);
%             outputState = applyToState(rev_cnot, outputState);
%             outputState = applyToState(cnot, outputState);
            sz = size(inputState);
            tmax = numQubits+1-min(obj.targetQubit, obj.targetQubit2);
            tmin = numQubits+1-max(obj.targetQubit, obj.targetQubit2);
            outputState = reshape(inputState, 2*ones(1, log2(numel(inputState))));
            outputState = permute(outputState, ...
                [1:tmin-1 tmax tmin+1:tmax-1 tmin tmax+1:ndims(outputState)]);
            outputState = reshape(outputState, sz);
        end

        function [instr, def, map] = getQASM(obj, ~, map, varargin)
            % Assumed Supported and never renamed
            instr = quantum.internal.gate.formatQASMInstruction(obj);
            def = "";
        end

        function qb = getQubits(obj)
            qb = [obj.targetQubit obj.targetQubit2];
        end

        function obj = setQubits(obj, qb)
            assert(length(qb) == 2);
            obj.targetQubit = qb(1);
            obj.targetQubit2 = qb(2);
        end

        function qb = getTargetQubits(obj)
            qb = [obj.targetQubit obj.targetQubit2];
        end

        function qb = getControlQubits(~)
            qb = zeros(1, 0);
        end
    end

    methods(Hidden, Static)
        function [Nctrl, Ntrgt, Nangle] = numInputs
            Nctrl = 0;
            Ntrgt = 2;
            Nangle = 0;
        end
    end
end
