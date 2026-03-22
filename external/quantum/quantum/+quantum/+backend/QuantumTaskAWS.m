classdef (Sealed) QuantumTaskAWS < quantum.backend.QuantumTask
    %QUANTUMTASKAWS A task sent to AWS for execution on a remote quantum device.
    %
    %   The quantumCircuit/run sends a task to run on a quantum device and
    %   returns a QuantumTaskAWS object. Additionally, a QuantumTaskAWS
    %   handle object can be attached to an existing remote task on AWS
    %   using that task's ARN (Amazon Resource Name):
    %
    %   TASK = QUANTUMTASKAWS(taskARN) returns a QuantumTaskAWS handle
    %   object that is attached to an existing remote task on AWS using
    %   that task's ARN. The task's ARN can be the previously saved TaskARN
    %   property of the QuantumTaskAWS object of a previously running task,
    %   or it can be retrieved from the AWS web interface.
    %
    %   Saving and loading of QuantumTaskAWS objects are not supported. To
    %   check on a task from a previous MATLAB session, save the
    %   QuantumTaskAWS object's TaskARN property from the previous session
    %   and construct a new QuantumTaskAWS object in the new MATLAB
    %   session.
    %
    %   Remote tasks expire after some time, we recommend saving the output
    %   of fetchOutput to keep the results available long-term.
    %
    %   Example:
    %       % Send a circuit to run on a remote quantum device
    %       % NOTE: This requires you to have access to AWS Braket, and
    %       % calling run sends a quantum task to AWS which will be charged
    %       % to your AWS account.
    %       gates = [hGate(1); cxGate(1, 2)];
    %       circ = quantumCircuit(gates);
    %       dev = quantum.backend.QuantumDeviceAWS("Lucy");
    %       task = run(circ, dev)
    %       task.Status
    %       wait(task)
    %       measurement = fetchOutput(task)
    %
    %   QuantumTaskAWS properties:
    %      Status           - Status of the task.
    %      TaskARN          - ARN device identifier.
    %
    %   QuantumTaskAWS methods:
    %      wait              - Wait for task to finish running or reach a specified status.
    %      fetchOutput       - Retrieve result of finished task.
    %      fetchDetails      - Retrieve details about the task.
    %      cancel            - Cancel the task.
    %
    %   See also quantum.backend.QuantumDeviceAWS, quantumCircuit/run

    %   Copyright 2022-2023 The MathWorks, Inc.


    properties (GetAccess = public, SetAccess = public)
        %TASKARN - The ARN identifier of the task
        %   TaskARN is the Amazon Resource Name (ARN) of the task.
        %
        %   This can be saved and used in a future MATLAB session to
        %   construct a new QuantumTaskAWS object and query the task's status.
        %
        %   See also quantum.backend.QuantumDeviceAWS
        TaskARN
    end

    properties(GetAccess=private, SetAccess=private)
        TaskRegion           % Extracted from TaskARN, stored separately for convenience.
        TaskDetailStruct     % getQuantumTask struct
        ResultDetailStr      % fileread S3 string
        WasCancelled = false % Set to true if task was cancelled
    end


    methods
        function obj = QuantumTaskAWS(taskARN)
            arguments
                taskARN {mustBeTextScalar}
            end

            if ismissing(taskARN)
                obj.TaskARN = taskARN;
                obj.Status = string(NaN);
                return
            end

            taskARN = string(taskARN);
            taskARNParts = split(taskARN, ':');
            if length(taskARNParts) ~= 6 || ~startsWith(taskARN, "arn:aws:braket") || ~startsWith(taskARNParts(6), "quantum-task")
                error(message("quantum:QuantumTaskAWS:InvalidTaskARN"))
            end

            obj.TaskARN = taskARN;
            obj.TaskRegion = taskARNParts(4);
            obj.updateTask();
        end

        function cancel(obj)
            %CANCEL  Cancel the task
            %   CANCEL(task) cancels the task on the AWS quantum computing
            %   service. If the task's Status is already "finished",
            %   cancelling has no effect.
            %
            %   Cancelling a task may cancel the cost of running the task
            %   depending on when it is cancelled, please consult the AWS
            %   Braket documentation for details.
            %
            %   See also quantum.backend.QuantumTaskAWS
            if ismissing(obj.TaskARN)
                error(message("quantum:QuantumTask:MissingNotSupported"))
            end

            obj.WasCancelled = true;

            status = obj.Status;
            if ismember(status, ["queued", "running"])
                quantum.internal.aws.cancelQuantumTask(obj.TaskARN, obj.TaskRegion);
                obj.updateTask();
            end
        end

        function s = fetchDetails(obj)
            %FETCHDETAILS Retrieve details about the task
            %
            %   s = FETCHDETAILS(task) returns a struct containing
            %   information about the task (in s.Task), and about its
            %   result if available (in s.Result). Each of these fields is
            %   again a struct containing various fields provided by the
            %   AWS quantum computing service.
            %
            %   See also quantum.backend.QuantumTaskAWS.

            arguments
                obj quantum.backend.QuantumTaskAWS
            end

            if ismissing(obj.TaskARN)
                error(message("quantum:QuantumTask:MissingNotSupported"))
            end

            s.Task = obj.TaskDetailStruct;

            status = obj.Status;
            if strcmp(status, "finished") && ~strcmp(obj.TaskDetailStruct.status, "CANCELLED")
                if isempty(obj.ResultDetailStr)
                    obj.ResultDetailStr = readS3Bucket(obj);
                end

                s.Result = jsondecode(obj.ResultDetailStr);
            end
        end

        function output = fetchOutput(obj)
            %FETCHOUTPUT Retrieve result of the finished task
            %   m = FETCHOUTPUT(task) returns a QuantumMeasurement object
            %   representing the output of the task.
            %
            %   If task is not in status "finished", fetchOutput will error.
            %
            %   See also quantum.backend.QuantumTaskAWS,
            %   quantum.gate.QuantumMeasurement.

            if ismissing(obj.TaskARN)
                error(message("quantum:QuantumTask:MissingNotSupported"))
            end

            status = obj.Status;
            switch status
                case {'queued', 'running'}
                    error(message("quantum:QuantumTaskAWS:OutputNotAvailable", status))
                case "failed"
                    error(message("quantum:QuantumTaskAWS:OutputNotAvailableFailed", obj.TaskDetailStruct.failureReason))
                case "finished"
                    if strcmp(obj.TaskDetailStruct.status, "CANCELLED")
                        error(message("quantum:QuantumTaskAWS:OutputNotAvailableCancelled"))
                    end

                    resultStr = readS3Bucket(obj);
                    results = decodeS3Contents(resultStr);
                    taskType = results.braketSchemaHeader.name;

                    if ~strcmp(taskType, "braket.task_result.gate_model_task_result")
                        error(message("quantum:QuantumTaskAWS:InvalidResultType"))
                    end

                    obj.ResultDetailStr = resultStr;
                    output = decodeGateResults(results);
            end
        end
    end

    methods(Access=protected)
        function updateTask(obj)
            if ismissing(obj.TaskARN)
                return
            end

            details = quantum.internal.aws.getQuantumTask(obj.TaskARN, obj.TaskRegion);
            
            obj.TaskDetailStruct = details;
            status = mapBraketStatusToMATLABStatus(details.status);
            obj.Status = status;

            if obj.WasCancelled && details.status == "QUEUED"
                % First cancellation seems to not have reached AWS, resend.
                quantum.internal.aws.cancelQuantumTask(obj.TaskARN, obj.TaskRegion);
            end
        end
    end

    methods (Hidden)
        function s = saveobj(~)
            s = struct;
            warning(message("quantum:QuantumTask:SaveNotSupported"));
        end
    end
    methods(Hidden, Static)
        function obj = loadobj(~)
            warning(message('quantum:QuantumTask:LoadNotSupported'));
            obj = quantum.backend.QuantumTaskAWS(string(NaN));
        end
    end
end



function content = decodeS3Contents(contentStr)
try
    content = jsondecode(contentStr);
catch causeME
    throw(addCause(MException(message("quantum:QuantumTaskAWS:InvalidS3Contents")), causeME))
end
end

function resultStr = readS3Bucket(obj)
s3bucket = obj.TaskDetailStruct.outputS3Bucket;
s3folder = obj.TaskDetailStruct.outputS3Directory;
jsonS3URI = "s3://"+s3bucket+"/"+s3folder+"/"+"results.json";
try
    resultStr = string(fileread(jsonS3URI));
catch causeME
    throw(addCause(MException(message("quantum:QuantumTaskAWS:FailedReadS3")), causeME))
end
end

function out = decodeGateResults(results)
% Cases are Braket result pragmas used in our QASM to standardize how the
% results are formatted.  

if isfield(results, "measurements")
    map = '01';
    [NumShots, NumQubits] = size(results.measurements);
    Obs = map(results.measurements+1);
    Obs = string(reshape(Obs, NumShots, NumQubits));
    [States, ~, ind] = unique(Obs);
    Counts = accumarray(ind, 1, size(States));
elseif isfield(results, "measurementProbabilities")
    % Shouldn't happen for tasks we send, keeping as a fall-back
    % for now
    States = string(fields(results.measurementProbabilities));
    States = erase(States, "x"); % caused by jsondecode numeric fieldnames
    Probs = struct2cell(results.measurementProbabilities);
    Probs = [Probs{:}]';
    NumShots = results.taskMetadata.shots;
    Counts = round(Probs*NumShots);
else
    error(message("quantum:QuantumTaskAWS:InvalidResultType"))
end
% Left-most bit corresponds to the top wire in the circuit.
out = quantum.gate.QuantumMeasurement(States, Counts);
end

function statusMATLAB = mapBraketStatusToMATLABStatus(statusBraket)
% CANCELLING is not a final status, might end in COMPLETED if task
% was already RUNNING
% CREATED might be returned by AWS after we create the task,
% so it doesnt match our PCT pending state.
switch statusBraket
    case  {'CANCELLING', 'QUEUED', 'CREATED'}
        statusMATLAB = "queued";
    case "RUNNING"
        statusMATLAB = "running";
    case {'COMPLETED','CANCELLED'}
        statusMATLAB = "finished";
    otherwise
        assert(statusBraket=="FAILED")
        statusMATLAB = "failed";
end
end
