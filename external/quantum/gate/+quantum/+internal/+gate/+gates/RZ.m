classdef RZ < quantum.internal.gate.gates.OneQubitGate
    %

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix
        theta
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "RZ";
        type = "rz";
    end

    methods
        function obj = RZ(targetQubit, theta)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
            obj.theta = theta;
            obj.gateMatrix =  [exp(-1i*theta/2), 0; 0, exp(1i*theta/2)];
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.RZ(obj.targetQubit, -obj.theta);
        end
    end

    methods (Hidden)
        function a = angles(obj)
            a = obj.theta;
        end

        function obj = setAngles(obj, theta)
            obj.theta = theta;
            obj.gateMatrix =  [exp(-1i*theta/2), 0; 0, exp(1i*theta/2)];
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
