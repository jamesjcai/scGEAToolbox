classdef Gate
    % This defines all methods that any gate needs to provide (most as abstract). The goal in the next steps is
    % for this to become the parent of all the Gate objects, instead of QubitGate being that parent.
    %
    % Currently some methods of quantumCircuit are just accessing internal properties of the gates. That should
    % be replaced with additional methods, which will be added here, too.

    %   Copyright 2021-2022 The MathWorks, Inc.

    properties(Abstract, SetAccess = private, Hidden, GetAccess = public)
        nameForPlot
    end

    methods(Abstract)
        % Apply the gate to a state vector
        outputState = applyToState(obj, inputState, numQubits)

        % Generate QASM code for this gate
        [instr, def, varargout] = getQASM(obj, varargin)

        qb = getQubits(obj)

        obj = setQubits(obj, qb)

        obj = inv(obj) % Could also call this ctranspose
    end

    methods(Hidden)

        function a = angles(~)
            a = zeros(1, 0);
        end

        function obj = setAngles(obj, ang)
            % No-op, any gate that has angles will overload this to
            % actually do something
            assert(isempty(ang));
        end
    end
end