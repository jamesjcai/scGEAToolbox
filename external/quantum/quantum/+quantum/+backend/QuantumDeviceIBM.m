classdef (Sealed) QuantumDeviceIBM < quantum.backend.QuantumDevice
    %QuantumDeviceIBM A quantum device available through IBM
    %
    %   DEV = QuantumDeviceIBM(deviceName) returns a representation of the
    %   specified device.
    %
    %   DEV = QuantumDeviceIBM(...,'FileName',fileName) specifies the
    %   fileName of the JSON file containing credentials used to
    %   authenticate with the Qiskit Runtime service. By default, fileName
    %   is "qiskit-ibm.json" in the .qiskit folder in the home directory.
    %
    %   DEV = QuantumDeviceIBM(...,'AccountName',accName) specifies the
    %   desired account to use from the JSON file. By default, the first
    %   account listed in fileName is used.
    %
    %   DEV = QuantumDeviceIBM(...,'UseSession',tf) toggles whether to
    %   create tasks in a session. A session is useful for reducing overall
    %   queue time when running multiple tasks, but provides no queue time
    %   benefit for a single task.
    %
    %   Saving and loading of QuantumDeviceIBM objects are not supported.
    %
    %   Example:
    %       % Construct a QuantumDeviceIBM object and run a circuit on it.
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
    %   QuantumDeviceIBM properties:
    %      Name             - Name of the device.
    %      AccountName      - Name of the account used for authorization.
    %      UseSession       - Toggles if a task should be created in a
    %                         session.
    %
    %   QuantumDeviceIBM methods:
    %      fetchDetails      - Retrieve details about the device
    %
    %   See also quantum.backend.QuantumTaskIBM, quantumCircuit/run

    %   Copyright 2022-2025 The MathWorks, Inc.
    properties(GetAccess=public, SetAccess=private)
        %NAME - Name of the device
        %
        %   See also quantum.backend.QuantumTaskIBM
        Name

        %ACCOUNTNAME - Name of the account used for authentication
        %
        %   The credentials listed under AccountName in the JSON file are
        %   used to connect to IBM. The cost of running a task is billed
        %   to this account.
        %
        %   See also quantum.backend.QuantumTaskIBM
        AccountName

        %USESESSION - Toggle for creating tasks in a session
        %
        %   A session is useful for reducing overall queue time when
        %   running multiple tasks, but provides no queue time benefit for
        %   a single task.
        %
        %   See also quantum.backend.QuantumTaskIBM
        UseSession
    end

    properties(Access = {?iQuantumDeviceIBM})
        SessionID = missing
    end

    properties(Access=private)
        Credentials
    end

    methods
        function obj = QuantumDeviceIBM(deviceName, NameValueArgs)
            arguments
                deviceName {mustBeTextScalar}
                NameValueArgs.UseSession {mustBeA(NameValueArgs.UseSession,'logical')} = false
                NameValueArgs.FileName {mustBeTextScalar} = quantum.internal.ibm.defaultFilepath("qiskit-ibm.json")
                NameValueArgs.AccountName {mustBeTextScalar} = ""
            end

            if ismissing(deviceName)
                obj.Name = string(NaN);
                obj.AccountName = string(NaN);
                obj.UseSession = false;
                return
            end

            deviceName = string(lower(deviceName));
            fileName = string(NameValueArgs.FileName);
            accountName = string(NameValueArgs.AccountName);

            try
                % Parse file for account name and the struct of credentials
                [credentials, accountName] = quantum.internal.ibm.readCredentials(fileName,accountName);
            catch ME
                throw(addCause(MException(message("quantum:QuantumDeviceIBM:FileReadError", fileName)), ME))
            end

            % Authenticate credentials by querying available devices
            credentials = quantum.internal.ibm.generateAccessToken(credentials);
            availableDevices = quantum.internal.ibm.listDevices(credentials);

            if ~ismember(deviceName, availableDevices)
                error(message("quantum:QuantumDeviceIBM:DeviceNotFound", accountName))
            end

            obj.Credentials = credentials;
            obj.Name = deviceName;
            obj.AccountName = accountName;

            if NameValueArgs.UseSession
                det = quantum.internal.ibm.getDevice(obj.Credentials, obj.Name);
                if det.Configuration.simulator
                    error(message("quantum:QuantumDeviceIBM:SessionNotSupported"))
                end
            end

            obj.UseSession = NameValueArgs.UseSession;
        end

        function details = fetchDetails(obj)
            %FETCHDETAILS Retrieve details about the quantum device
            %
            %   s = FETCHDETAILS(dev) returns a struct containing
            %   information about the device.
            %
            %   The output of this function may change over time. Calling
            %   FETCHDETAILS always retrieves the latest version provided
            %   by IBM.
            %
            %   See also quantum.backend.QuantumDeviceIBM.

            if ismissing(obj.Name)
                error(message("quantum:QuantumDevice:MissingNotSupported"))
            end

            details = getDetails(obj);
        end
    end

    methods(Hidden)
        function task = sendCircuit(obj, circuits, NameValueArgs)
            arguments
                obj quantum.backend.QuantumDeviceIBM
                circuits quantumCircuit
                NameValueArgs.NumShots (1,1) {mustBeInteger, mustBePositive} = 100
                % NOTE: OptimizationLevel was deprecated by IBM. It is no longer 
                % sent to the server.
                NameValueArgs.UseErrorMitigation (1,1) {mustBeA(NameValueArgs.UseErrorMitigation, 'logical')} = true
                NameValueArgs.Observable
            end
            
            [qasms, finalMaps] = transpileQASM(obj, circuits);

            if isfield(NameValueArgs, "Observable")
                obs = NameValueArgs.Observable;
                if ~isa(obs, "observable")
                    error(message('quantum:observable:invalidObservable'));
                end
                % Check expansion
                nC = [circuits.NumQubits];
                nO = [obs.NumQubits];
                [cIdx,oIdx] = quantum.internal.gate.setupLocalImplicitExpansion(size(qasms), size(obs));
                nC = nC(cIdx);
                nO = nO(oIdx);
                if ~isequal(nC(:), nO(:))
                    error(message("quantum:observable:observeIncorrectNumQubits"))
                end
                NameValueArgs.FinalMaps = finalMaps;
            end

            task = sendQASM(obj, qasms, NameValueArgs);
        end

        function [qasms, finalMaps] = transpileQASM(obj, circuits)

            details = fetchDetails(obj);

            if details.Configuration.simulator
                % Transpilation must be used for simulator to decompose
                % gates but they can be applied between any qubits.
                adj = true(details.Configuration.n_qubits);
                basisGates = string(details.Configuration.basis_gates);
                if all(ismember(["rz" "sx" "x" "cz"], basisGates))
                    % Simulator is expected to always have this basis set.
                    mode = 6;
                else
                    error(message("quantum:generateQASM:missingSupportedGates"))
                end
                % The map is not added as code comment here because error
                % mitigation is only used for non-simulators.
                [qasms, finalMaps] = quantum.internal.transpile.transpile(circuits, adj, mode);
                return
            end

            S = details.Configuration.gates;
            T = struct2table(S);
            if ~all(isfield(S, ["coupling_map" "name"]))
                error(message("quantum:generateQASM:missingSupportedGates"))
            end

            basisGates = lower(string(T.name));

            if all(ismember(["rz" "sx" "x" "ecr"], basisGates))
                % Eagle processor
                mode = 5;
                idx = matches(basisGates, "ecr");
            elseif all(ismember(["rz" "sx" "x" "cz"], basisGates))
                % Heron processor
                mode = 6;
                idx = matches(basisGates, "cz");
            else
                error(message("quantum:generateQASM:missingSupportedGates"))
            end

            % If the two-qubit basis gate is bidirectional, both edge
            % directions are returned from the server.
            supportedEdges = T.coupling_map{idx}+1;
            G = digraph(supportedEdges(:, 1), supportedEdges(:, 2));
            adj = logical(full(adjacency(G)));

            [qasms, finalMaps] = quantum.internal.transpile.transpile(circuits, adj, mode);

            if mode==5
                % Add empty ECR definition to all qasm codes to pass server
                % checks
                ecrDef = newline+"gate ecr q0,q1 {}";
                qasms = insertAfter(qasms, 'include "stdgates.inc";', ecrDef);
            end

            L = numel(qasms);
            for ii = 1:L
                % Add each final map as code comment. This is parsed by the
                % task when fetching results for error mitigation.
                qasms(ii) = insertVectorCommentAfterHeader(qasms(ii), finalMaps{ii}, "finalMap");
            end
        end

        function task = sendQASM(obj, qasms, NameValueArgs)
            if ismissing(obj.Name)
                error(message('quantum:QuantumDevice:MissingNotSupported'));
            end
            
            taskInfo.backend = obj.Name;
            runInfo = struct;
            
            % The same NumShots is used for each circuit
            if isfield(NameValueArgs, "Observable")
                taskInfo.program_id = "estimator";
                pubs = setupEstimatorPUB(qasms, NameValueArgs);
                % NumShots is added outside the PUB.
                runInfo.options.default_shots = NameValueArgs.NumShots;
            else
                taskInfo.program_id = "sampler";
                % NumShots is added within the PUB
                pubs = setupSamplerPUB(qasms, NameValueArgs);
            end

            runInfo.pubs = pubs;
            
            runInfo.version = 2;
            runInfo.support_qiskit = false;
            
            if NameValueArgs.UseErrorMitigation
                % Enable error mitigation features when running the circuit
                runInfo.options.twirling.enable_gates = true;
                runInfo.options.twirling.enable_measure = true;
                % Enable post processing to mitigate measurement errors.
                % This is done by the task when fetching the output.
                taskInfo.tags = {'{"useErrorMitigation":true,"version":2}'};
            else
                taskInfo.tags = {'{"useErrorMitigation":false,"version":2}'};
            end

            taskInfo.params = runInfo;

            if ~obj.UseSession
                % Create task without session
                task = send(obj, taskInfo);
                return
            end

            % Struct of session options
            sessionInfo.backend = obj.Name;
            sessionInfo.mode = "batch";

            if ismissing(obj.SessionID)
                % Create new session
                obj.SessionID = quantum.internal.ibm.createSession(obj.Credentials, sessionInfo);
            else
                % Use the current session if it's open otherwise create
                % new session
                sessionDetails = quantum.internal.ibm.getSession(obj.Credentials, obj.SessionID);
                isSessionOpen = sessionDetails.accepting_jobs && ~strcmpi(sessionDetails.state, 'closed');
                if ~isSessionOpen
                    % Create new session. The old session will be closed by the
                    % server when all tasks are finished
                    obj.SessionID = quantum.internal.ibm.createSession(obj.Credentials, sessionInfo);
                end
            end

            % Create task with session
            taskInfo.session_id = obj.SessionID;
            task = send(obj, taskInfo);
        end
    end

    methods(Hidden)
        function s = saveobj(~)
            s = struct;
            warning(message("quantum:QuantumDevice:SaveNotSupported"));
        end
    end

    methods (Hidden, Static)
        function obj = loadobj(~)
            warning(message('quantum:QuantumDevice:LoadNotSupported'));
            obj = quantum.backend.QuantumDeviceIBM(string(NaN));
        end
    end

    methods(Hidden, Access = {?iQuantumDeviceIBM})
        function closeSession(obj)
            quantum.internal.ibm.closeSession(obj.Credentials, obj.SessionID);
        end
    end

    methods(Access=private)
        function details = getDetails(obj)
            details = quantum.internal.ibm.getDevice(obj.Credentials, obj.Name);
        end

        function task = send(obj, taskInfo)
            taskID = quantum.internal.ibm.createQuantumTask(obj.Credentials, taskInfo);
            task = quantum.backend.QuantumTaskIBM.createFromCredentials(taskID, obj.AccountName, obj.Credentials);
        end
    end
end

% The IBM Runtime interface accepts Primitive-Unified-Blocks (PUBs) that
% contains circuits and additional info. Each 1-by-3 cell below represents
% one PUB. The parameters are left empty because our circuits are 
% non-symbolic.

function pubs = setupSamplerPUB(qasms, NameValueArgs)
% Sampler Interface PUB:
% {circuit, parameters, numShots}

szOut = size(qasms);

% Add permutation for this group of circuits to the top code.
% This is parsed by the task when fetching results to return the expected size.
qasms(1) = insertVectorCommentAfterHeader(qasms(1), szOut, "szOut");

d1 = numel(qasms);
pubs = cell(d1,1);
for ii = 1:d1
    % This PUB contains 1 circuit
    pubs{ii} = {qasms(ii), [], NameValueArgs.NumShots};
end
end


function pubs = setupEstimatorPUB(qasms, NameValueArgs)
% Estimator Interface PUB:
% {circuit, observables, parameters}

obs = NameValueArgs.Observable;
finalMaps = NameValueArgs.FinalMaps;

[oIdx, perm, szOutIdxPerm] = quantum.internal.gate.setupRemoteImplicitExpansion(size(qasms), size(obs));

% Add permutations for this group of circuits/observables to the top code.
% This is parsed by the task when fetching results to return the expected size.
qasms(1) = insertVectorCommentAfterHeader(qasms(1), perm, "perm");
qasms(1) = insertVectorCommentAfterHeader(qasms(1), szOutIdxPerm, "szOutIdxPerm");

% Implicit expansion
d1 = numel(qasms);
d2 = size(oIdx,2);
pubs = cell(d1,1);
for ii=1:d1

    S = cell(1,d2);
    for jj=1:d2

        o = obs(oIdx(ii,jj));
        p = char(o.Paulis);
        w = o.Weights;

        % Permute with final map and reverse to match convention
        map = finalMaps{ii}+1;
        mp = repmat('I', [length(w) max(map)]);
        mp(:, map) = fliplr(p);

        % Build struct for JSON
        s = struct;
        for kk = 1:size(mp,1)
            s.(mp(kk,:)) = w(kk);
        end
        S{jj} = s;
    end
        
    % This PUB contains 1 circuit with the expanded observables
    pubs{ii} = {qasms(ii), S, {}};
end
end


function qasm = insertVectorCommentAfterHeader(qasm, vec, name)
header = 'include "stdgates.inc";';
str = join(string(vec), ',');
comment = newline+"//Start "+name+newline+...
    "//"+str+newline+...
    "//End "+name;
qasm = insertAfter(qasm, header, comment);
end