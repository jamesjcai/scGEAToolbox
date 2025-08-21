classdef ControlledQubitGate < quantum.internal.gate.gates.Gate
    %

    %   Copyright 2021-2022 The MathWorks, Inc.

    properties(Abstract, SetAccess = protected, GetAccess = public)
        controlQubit
        targetQubit
        gateMatrix
    end

    properties(Abstract, SetAccess = private, Hidden, GetAccess = public)
        type
    end

    methods
        %constructor
        function obj = ControlledQubitGate(controlQubit, targetQubit)
            if targetQubit == controlQubit
                error(message('quantum:gates:matchingQubits'));
            end
            obj.targetQubit = targetQubit;
            obj.controlQubit = controlQubit;
        end

        function outputState = applyToState(obj, inputState, numQubits)
            if numel(inputState) < 4
                error(message('quantum:gates:tooFewQubits', 2));
            end
            state0 = [1 0; 0 0];
            state1 = [0 0; 0 1];

            temp0 = quantum.internal.gate.applyMatToDim(inputState, state0, numQubits+1-obj.controlQubit);
            temp1 = quantum.internal.gate.applyMatToDim(inputState, state1, numQubits+1-obj.controlQubit);
            temp1 = quantum.internal.gate.applyMatToDim(temp1, obj.gateMatrix, numQubits+1-obj.targetQubit);

            outputState = temp0 + temp1;

%             y1 = g2.applyToState(x);
%             y = x;
%             y(2, :) = y1(2, :);
        end
        
        function qb = getQubits(obj)
            qb = [obj.controlQubit obj.targetQubit];
        end

        function obj = setQubits(obj, qb)
            assert(length(qb) == 2);
            obj.controlQubit = qb(1);
            obj.targetQubit = qb(2);
        end

        function qb = getTargetQubits(obj)
            qb = obj.targetQubit;
        end

        function qb = getControlQubits(obj)
            qb = obj.controlQubit;
        end
    end
end
