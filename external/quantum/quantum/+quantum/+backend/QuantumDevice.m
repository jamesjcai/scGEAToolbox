classdef QuantumDevice < handle & matlab.mixin.Scalar
    % Any child of this class is accepted by quantumCircuit/run, which will
    % call its hidden method:
    %
    %    task = sendCircuit(device, circuit, options)
    %
    % The return value of this must be a class that inherits from
    % quantum.backend.QuantumTask.
    %
    % See also quantum.backend.QuantumDeviceAWS

    %   Copyright 2023 The MathWorks, Inc.

    methods (Abstract)
        s = fetchDetails(device)
    end

    methods (Abstract, Hidden)
        task = sendCircuit(device, circuit, options)
    end
end