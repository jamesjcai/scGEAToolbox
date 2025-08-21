classdef Missing < quantum.internal.gate.gates.Gate
    % This is used to represent a missing gate. These will arise from
    % indexing expressions that need a QuantumGate array to be padded with
    % default elements, or for a gate from a future version which doesn't
    % exist in this version - when loading such a gate from a .mat file, we
    % represent it by this.

    %   Copyright 2021-2022 The MathWorks, Inc.


    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "Missing";
        type = "missing";
    end

    methods
        %constructor
        function obj = Missing()
            
        end

        function obj = inv(~) %#ok<STOUT>
            error(message('quantum:gates:MissingGateNotSupported'))
        end
    end

    methods(Hidden)
        function outputState = applyToState(~, ~, ~) %#ok<STOUT>
            error(message('quantum:gates:MissingGateNotSupported'))
        end

        function [instr, def, map] = getQASM(~, varargin) %#ok<STOUT>
            error(message('quantum:gates:MissingGateNotSupported'))
        end

        function qb = getQubits(~)
            qb = zeros(1, 0);
        end

        function obj = setQubits(~, ~) %#ok<STOUT>
            error(message('quantum:gates:MissingGateNotSupported'))
        end

        function qb = getTargetQubits(~)
            qb = zeros(1, 0);
        end

        function qb = getControlQubits(~)
            qb = zeros(1, 0);
        end
    end
end
