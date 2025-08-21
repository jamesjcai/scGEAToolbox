classdef RY < quantum.internal.gate.gates.OneQubitGate
    %

    %   Copyright 2021-2023 The MathWorks, Inc.
    
    properties(SetAccess = protected, GetAccess = public)
        targetQubit
        gateMatrix
        theta
    end

    properties(SetAccess = private, Hidden, GetAccess = public)
        nameForPlot = "RY";
        type = "ry";
    end

    methods
        function obj = RY(targetQubit, theta)
            obj@quantum.internal.gate.gates.OneQubitGate(targetQubit);
            obj.theta = theta;
            obj.gateMatrix = [cos(theta/2), -sin(theta/2);...
                sin(theta/2), cos(theta/2)];
        end

        function obj = inv(obj)
            obj = quantum.internal.gate.gates.RY(obj.targetQubit, -obj.theta);
        end
    end

    methods (Hidden)
        function a = angles(obj)
            a = obj.theta;
        end

        function obj = setAngles(obj, theta)
            obj.theta = theta;
            obj.gateMatrix = [cos(theta/2), -sin(theta/2);...
                sin(theta/2), cos(theta/2)];
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
