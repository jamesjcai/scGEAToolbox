classdef ControlledControlledQubitGate < quantum.internal.gate.gates.Gate
    %

    %   Copyright 2021-2022 The MathWorks, Inc.

    properties(Abstract, SetAccess = protected, GetAccess = public)
        controlQubit1
        controlQubit2
        targetQubit
        gateMatrix
    end

    properties(Abstract, SetAccess = private, Hidden, GetAccess = public)
        type
    end

    methods
        %constructor
        function obj = ControlledControlledQubitGate(controlQubit1,...
                controlQubit2, targetQubit)
            if targetQubit == controlQubit1 || targetQubit == controlQubit2 || ...
                    controlQubit1 == controlQubit2
                error(message('quantum:gates:matchingQubits'));
            end
            obj.targetQubit = targetQubit;
            obj.controlQubit1 = controlQubit1;
            obj.controlQubit2 = controlQubit2;
        end

        function outputState = applyToState(obj, inputState, numQubits)
            if numel(inputState) < 8
                error(message('quantum:gates:tooFewQubits', 3));
            end
            state1 = [0 0; 0 1];

            controlledState = inputState;
            controlledState = quantum.internal.gate.applyMatToDim(controlledState, state1, ...
                numQubits+1-obj.controlQubit1);
            controlledState = quantum.internal.gate.applyMatToDim(controlledState, state1, ...
                numQubits+1-obj.controlQubit2);

            outputState = inputState - controlledState + ...
                quantum.internal.gate.applyMatToDim(controlledState, obj.gateMatrix, ...
                numQubits+1-obj.targetQubit);
        end

        function qb = getQubits(obj)
            qb = [obj.controlQubit1 obj.controlQubit2 obj.targetQubit];
        end

        function obj = setQubits(obj, qb)
            assert(length(qb) == 3);
            obj.controlQubit1 = qb(1);
            obj.controlQubit2 = qb(2);
            obj.targetQubit = qb(3);
        end

        function qb = getTargetQubits(obj)
            qb = obj.targetQubit;
        end

        function qb = getControlQubits(obj)
            qb = [obj.controlQubit1 obj.controlQubit2];
        end
    end
end
