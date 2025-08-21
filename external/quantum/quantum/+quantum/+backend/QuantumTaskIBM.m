classdef (Sealed) QuantumTaskIBM < quantum.backend.QuantumTask
    %QUANTUMTASKIBM A task sent to IBM for execution on a remote quantum device.
    %
    %   The quantumCircuit/run method sends a task to run on a quantum
    %   device and returns a QuantumTaskIBM object. Additionally, a
    %   QuantumTaskIBM handle object can be created from an existing remote
    %   task on IBM using that task's ID.
    %
    %   TASK = QuantumTaskIBM(taskID) returns a QuantumTaskIBM handle
    %   object that is created from an existing remote task on IBM using
    %   that task's ID. The task's ID can be the previously saved TaskID
    %   property of the QuantumTaskIBM object of a previously running task,
    %   or it can be retrieved from the IBM web interface.

    %   TASK = QuantumTaskIBM(...,'FileName',fileName) specifies the
    %   fileName of the JSON file containing credentials used to
    %   authenticate with the Qiskit Runtime service. By default, fileName
    %   is "qiskit-ibm.json" in the .qiskit folder in the home directory.
    %
    %   TASK = QuantumTaskIBM(...,'AccountName',accName) specifies the
    %   desired account to use from the JSON file. By default, the first
    %   account listed in fileName is used.
    %
    %   Saving and loading of QuantumTaskIBM objects are not supported. To
    %   check on a task from a previous MATLAB session, save the
    %   QuantumTaskIBM object's TaskID property from the previous session
    %   and construct a new QuantumTaskIBM object in the new MATLAB
    %   session.
    %
    %   Remote tasks expire after some time. Save the output of fetchOutput
    %   to keep the results available long-term.
    %
    %   Example:
    %       % Send a circuit to run on a remote quantum device.
    %       % NOTE: This requires you to have access to Qiskit Runtime, and
    %       % calling run sends a quantum task to IBM which will be charged
    %       % to your account.
    %       gates = [hGate(1); cxGate(1, 2)];
    %       circ = quantumCircuit(gates);
    %       dev = quantum.backend.QuantumDeviceIBM("ibmq_belem");
    %       task = run(circ, dev)
    %       task.Status
    %       wait(task)
    %       measurement = fetchOutput(task)
    %
    %   QuantumTaskIBM properties:
    %      Status          - Status of the task.
    %      TaskID          - Task identifier.
    %      SessionID       - Session identifier.
    %      AccountName     - Name of the account used for authorization.
    %
    %   QuantumTaskIBM methods:
    %      wait              - Wait for task to finish running or reach a specified status.
    %      fetchOutput       - Retrieve result of finished task.
    %      fetchDetails      - Retrieve details about the task.
    %      cancel            - Cancel the task.
    %
    %   See also quantum.backend.QuantumDeviceIBM, quantumCircuit/run

    %   Copyright 2022-2024 The MathWorks, Inc.
    properties (GetAccess=public, SetAccess=private)
        %TASKID - The identifier of the task
        %
        %   This can be saved and used in a future MATLAB session to
        %   construct a new QuantumTaskIBM object and query the task's status.
        %
        %   See also quantum.backend.QuantumDeviceIBM
        TaskID

        %SESSIONID- The identifier of the session
        %
        %   Defaults to <missing> when the task does not belong to a
        %   session.
        %
        %   See also quantum.backend.QuantumDeviceIBM
        SessionID = missing

        %ACCOUNTNAME - Name of the account used for authentication
        %
        %   The credentials listed under AccountName in the JSON file are
        %   used to connect to IBM. The cost of running a task is billed
        %   to this account.
        %
        %   See also quantum.backend.QuantumDeviceIBM
        AccountName
    end

    properties(Access=private)
        TaskDetails
        ResultDetails

        Credentials
    end

    methods
        function obj = QuantumTaskIBM(taskID, NameValueArgs)
            arguments
                taskID {mustBeTextScalar}
                NameValueArgs.FileName {mustBeTextScalar} = quantum.internal.ibm.defaultFilepath("qiskit-ibm.json")
                NameValueArgs.AccountName {mustBeTextScalar} = ""
            end

            if ismissing(taskID)
                obj.TaskID = string(NaN);
                obj.SessionID = string(NaN);
                obj.AccountName = string(NaN);
                obj.Status = string(NaN);
                return
            end

            taskID = string(taskID);
            fileName = string(NameValueArgs.FileName);
            accountName = string(NameValueArgs.AccountName);

            try
                % Parse file for account name and the struct of credentials
                [credentials, accountName] = quantum.internal.ibm.readCredentials(fileName,accountName);
            catch ME
                throw(addCause(MException(message("quantum:QuantumTaskIBM:FileReadError", fileName)), ME))
            end

            try
                % Authenticate credentials by querying task details
                credentials = quantum.internal.ibm.generateAccessToken(credentials);
                details = quantum.internal.ibm.getQuantumTaskDetails(credentials, taskID);
            catch ME
                serverMsg = ME.cause{1}.message;
                if contains(serverMsg, "404")
                    % Credentials are valid but the taskID was not found
                    % on the account
                    error(message("quantum:QuantumTaskIBM:TaskNotFound", accountName))
                else
                    throw(ME)
                end
            end

            obj.Credentials = credentials;

            obj.TaskID = string(details.id);
            obj.AccountName = accountName;
            obj.updateTask();
        end

        function cancel(obj)
            %CANCEL  Cancel the task
            %   CANCEL(task) cancels the task on the IBM quantum computing
            %   service. If the task's Status is already "finished",
            %   cancelling has no effect.
            %
            %   Cancelling a task may cancel the cost of running the task
            %   depending on when it is cancelled, please consult the Qiskit
            %   Runtime documentation for details.
            %
            %   See also quantum.backend.QuantumTaskIBM

            if ismissing(obj.TaskID)
                error(message("quantum:QuantumTask:MissingNotSupported"))
            end

            status = obj.Status;
            if ismember(status, ["queued", "running"])
                quantum.internal.ibm.cancelQuantumTask(obj.Credentials, obj.TaskID)
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
            %   Qiskit Runtime service.
            %
            %   See also quantum.backend.QuantumTaskIBM.

            if ismissing(obj.TaskID)
                error(message("quantum:QuantumTask:MissingNotSupported"))
            end

            % Query the status to get the latest task details
            status = obj.Status;
            s.Task = obj.TaskDetails;

            if isequal(status, "finished") && ~strcmpi(obj.TaskDetails.status, "Cancelled")
                if isempty(obj.ResultDetails)
                    obj.ResultDetails = quantum.internal.ibm.getQuantumTaskResults(obj.Credentials, obj.TaskID);
                end

                s.Result = obj.ResultDetails;
            end
        end

        function meas = fetchOutput(obj)
            %FETCHOUTPUT Retrieve result of the finished task
            %   meas = FETCHOUTPUT(task) returns a QuantumMeasurement object
            %   representing the output of the task.
            %
            %   If task is not in status "finished", fetchOutput will error.
            %
            %   See also quantum.backend.QuantumTaskIBM,
            %   quantum.gate.QuantumMeasurement.

            if ismissing(obj.TaskID)
                error(message("quantum:QuantumTask:MissingNotSupported"))
            end

            status = obj.Status;
            switch status
                case {'queued', 'running'}
                    error(message("quantum:QuantumTaskIBM:OutputNotAvailable", status))
                case "failed"
                    error(message("quantum:QuantumTaskIBM:OutputNotAvailableFailed", obj.TaskDetails.state.reason))
                case "finished"
                    if strcmpi(obj.TaskDetails.status, "Cancelled")
                        error(message("quantum:QuantumTaskIBM:OutputNotAvailableCancelled"))
                    end

                    if ~strcmpi(obj.TaskDetails.program.id, "sampler")
                        error(message("quantum:QuantumTaskIBM:InvalidResultType"))
                    end

                    results = quantum.internal.ibm.getQuantumTaskResults(obj.Credentials, obj.TaskID);

                    [meas, version, useErrorMitigation] = decodeGateResults(results, obj.TaskDetails);

                    obj.ResultDetails = results;

                    deviceName = string(obj.TaskDetails.backend);
                    det = quantum.internal.ibm.getDevice(obj.Credentials, deviceName);

                    if version==2 && useErrorMitigation && ~det.Configuration.simulator
                        % Extract device qubits used to run circuit.
                        try
                            % Task object expects pubs to be a scalar cell
                            % containing 1x3 cell {'OpenQASM 3.0; ...', [], numShots}
                            pubs = obj.TaskDetails.params.pubs;
                            qasm = string(pubs{1}{1});
                            mapStr = extractBetween(qasm, "//Start finalMap", "//End finalMap");
                            finalMap = str2num(erase(mapStr, "//")) + 1; %#ok<ST2NM>
                        catch ME
                            % Server format is different than expected
                            error(message("quantum:QuantumTaskIBM:InvalidResultType"))
                        end

                        % Get calibration matrices for these qubits
                        allCals = quantum.internal.ibm.getDeviceCalibration(obj.Credentials, deviceName);
                        cals = allCals(:,:,finalMap);

                        % Solve for mitigated probability
                        states = meas.MeasuredStates;
                        A = quantum.internal.gate.assignmentMatrix(states, cals);
                        mprobs = A\meas.Probabilities;
                        meas = quantum.gate.QuantumMeasurement(states, mprobs, "probs");
                    end
            end
        end
    end

    methods(Access=protected)
        function updateTask(obj)

            if ismissing(obj.TaskID)
                return
            end

            details = quantum.internal.ibm.getQuantumTaskDetails(obj.Credentials, obj.TaskID);

            obj.TaskDetails = details;
            obj.Status = mapIBMStatusToMATLABStatus(details.status);

            if isfield(details, "session_id") && ~isempty(details.session_id)
                % Only update the SessionID property when the task was ran
                % in a session
                obj.SessionID = string(details.session_id);
            end
        end
    end

    methods(Hidden)
        function s = saveobj(~)
            s = struct;
            warning(message("quantum:QuantumTask:SaveNotSupported"));
        end
    end

    methods(Hidden, Static)
        function obj = loadobj(~)
            warning(message('quantum:QuantumTask:LoadNotSupported'));
            obj = quantum.backend.QuantumTaskIBM(string(NaN));
        end

        function obj = createFromCredentials(taskID, accountName, credentials)
            obj = quantum.backend.QuantumTaskIBM(string(nan));
            obj.TaskID = taskID;
            obj.AccountName = accountName;
            obj.Credentials = credentials;
            obj.updateTask()
        end
    end
end

function [meas, version, useErrorMitigation] = decodeGateResults(results, details)
% Tasks may contain muliple results when created from another source. This
% object represents 1 task so only return first result.
resultIdx = 1;
if ismember('quasi_dists', fieldnames(results))
    % V1 Primitive returns error mitigated probabilities
    version = 1;
    useErrorMitigation = false;

    data = results.quasi_dists(resultIdx);
    states = reverse(string(erase(fields(data), "x")));
    probs = cell2mat(struct2cell(data));
    meas = quantum.gate.QuantumMeasurement(states, probs, "probs");
elseif ismember('results', fieldnames(results))
    % V2 Primitive returns raw bitstring samples
    tags = jsondecode(details.tags{:});
    version = tags.version;
    useErrorMitigation = tags.useErrorMitigation;

    data = results.results(resultIdx);
    samplesHex = string(data.data.c.samples);
    [statesHex, ~, ind] = unique(samplesHex);
    counts = accumarray(ind, 1, size(statesHex));
    numQubits = data.data.c.num_bits;
    states = quantum.internal.ibm.hex2bin(statesHex, numQubits);
    states = reverse(states);
    meas = quantum.gate.QuantumMeasurement(states, counts);
else
    error(message("quantum:QuantumTaskIBM:InvalidResultType"))
end
end

function statusMATLAB = mapIBMStatusToMATLABStatus(statusIBM)
switch statusIBM
    case "Queued"
        statusMATLAB = "queued";
    case "Running"
        statusMATLAB = "running";
    case {'Completed','Cancelled', 'Cancelled - Ran too long'}
        statusMATLAB = "finished";
    otherwise
        assert(statusIBM=="Failed")
        statusMATLAB = "failed";
end
end
