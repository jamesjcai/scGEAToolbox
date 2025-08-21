classdef QuantumTask < handle & matlab.mixin.Scalar
    % A child of this class is returned by quantumCircuit/run. It provides
    % access to a quantum task which is running on an external device.
    %
    % See also quantum.backend.QuantumTaskAWS

    %   Copyright 2023 The MathWorks, Inc.

    properties (GetAccess = public, SetAccess = protected)
        %STATUS - Status of the task
        %   Status is a string representing the current status of the task.
        %   Status is one of "queued", "running", "finished", or "failed".
        %
        %   Every check of the Status property causes a remote query to
        %   update it.
        %
        %   See also quantum.backend.QuantumTaskAWS
        Status
    end

    methods
        function status = get.Status(obj)
            obj.updateTask(); % retrieve latest information
            status = obj.Status;
        end
    end

    methods (Abstract)
        out = fetchOutput(task)
        s = fetchDetails(task)
        cancel(task)
    end

    methods (Abstract, Access = protected)
        updateTask
    end

    methods
        function tfout = wait(obj, status, timeout)
            %WAIT  Wait for task to reach final status
            %   WAIT(task) waits for the task to be "finished" (or
            %   "failed"). Internally, this method checks the Status
            %   property in a loop.
            %
            %   WAIT(task, status) waits for the task to reach a specified
            %   status. This status must be either "finished" or "running".
            %   If a task reaches the status "failed", wait stops waiting.
            %   The default status to wait for is "finished".
            %
            %   tf = WAIT(task, status, timeout) waits for the task to
            %   reach a specified status, for at most timeout seconds. The
            %   output tf is false if WAIT timed out early. Otherwise, tf
            %   is true. The default timeout is Inf.
            %
            %   See also quantum.backend.QuantumTaskAWS

            arguments
                obj quantum.backend.QuantumTask
                status {mustBeMember(status, ["finished", "running"])} = "finished"
                timeout (1,1) {mustBePositive, mustBeNumeric} = Inf
            end

            if ismissing(obj.Status)
                error(message("quantum:QuantumTask:MissingNotSupported"))
            end

            % allows for comparison; copies parallel.internal.types.States
            StatusEnum = dictionary(["queued","running","finished","failed"],[2,3,4,101]);
            startTime = tic;
            while true
                tf = StatusEnum(obj.Status) >= StatusEnum(status);
                if tf || toc(startTime) > timeout
                    break
                end
                pause(0.01);
            end

            if nargout > 0
                tfout = tf;
            end
        end
    end

end