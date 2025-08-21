classdef TwoQubitRotationGate < quantum.internal.gate.gates.Gate

    %   Copyright 2022 The MathWorks, Inc.

    properties(Abstract, SetAccess = protected, GetAccess = public)
        targetQubit1
        targetQubit2
        gateMatrix
    end

    methods
        %constructor
        function obj = TwoQubitRotationGate(targetQubit1, targetQubit2)
            if targetQubit1 == targetQubit2
                error(message('quantum:gates:matchingQubits'));
            end
            obj.targetQubit1 = targetQubit1;
            obj.targetQubit2 = targetQubit2;
        end

        function outputState = applyToState(obj, inputState, numQubits)
            if numel(inputState) < 4
                error(message('quantum:gates:tooFewQubits', 2));
            end

            sz = size(inputState);
            q1 = numQubits+1-obj.targetQubit1;
            q2 = numQubits+1-obj.targetQubit2;
            mind = max(q1, q2);
            innerSize = [2*ones(1, mind) numel(inputState)/2^mind];
            inputState = reshape(inputState, innerSize);
            perm = [q1 q2 setdiff(1:ndims(inputState), [q1 q2])];

            inputState = permute(inputState, perm);
            inputState = obj.gateMatrix * reshape(inputState, 4, []);
            inputState = reshape(inputState, innerSize);
            inputState = ipermute(inputState, perm);
            outputState = reshape(inputState, sz);
        end

        function qb = getQubits(obj)
            qb = [obj.targetQubit1 obj.targetQubit2];
        end

        function obj = setQubits(obj, qb)
            assert(length(qb) == 2);
            obj.targetQubit1 = qb(1);
            obj.targetQubit2 = qb(2);
        end

        function qb = getTargetQubits(obj)
            qb = [obj.targetQubit1 obj.targetQubit2];
        end

        function qb = getControlQubits(obj)
            qb = zeros(1, 0);
        end
    end
end
