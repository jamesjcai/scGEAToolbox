classdef (Sealed) QuantumDeviceAWS < quantum.backend.QuantumDevice
    %QuantumDeviceAWS A quantum device available through AWS
    %
    %   DEV = QuantumDeviceAWS(deviceName) returns a representation of the
    %   specified device, if it is uniquely defined by deviceName.
    %
    %   DEV = QuantumDeviceAWS(deviceARN) returns a representation of the
    %   device based on its ARN identifier. If deviceARN exists in multiple
    %   regions, the AWS_DEFAULT_REGION environment variable is used as the
    %   region.
    %
    %   DEV = QuantumDeviceAWS(...,'Region',region) specifies the region of
    %   the device. This can be used to uniquely identify a device if its
    %   name or ARN exists in multiple regions.
    %
    %   DEV = QuantumDeviceAWS(...,'S3Path',p) specifies the path in an S3
    %   bucket where the results will be written to and retrieved from. By
    %   default, S3Path is "s3://amazon-braket-mathworks/default". The
    %   specified path must start with "s3://amazon-braket-" and must
    %   define at least one folder inside of the bucket. This bucket must
    %   be accessible from your AWS account.
    %
    %   Saving and loading of QuantumDeviceAWS objects are not supported.
    %
    %   Example:
    %       % Construct a QuantumDeviceAWS object and run a circuit on it.
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
    %   QuantumDeviceAWS properties:
    %      Name             - Short device name.
    %      DeviceARN        - ARN device identifier.
    %      Region           - Region the device is in.
    %      S3Path           - Remote location where results are stored.
    %
    %   QuantumDeviceAWS methods:
    %      fetchDetails      - Retrieve details about the device
    %
    %   See also quantum.backend.QuantumTaskAWS, quantumCircuit/run

    %   Copyright 2022-2025 The MathWorks, Inc.

    properties (GetAccess = public, SetAccess = private)
        %NAME - Name of the device
        %   Name is a string, giving a short name for the quantum device.
        %
        %   See also quantum.backend.QuantumDeviceAWS
        Name

        %DEVICEARN - The ARN identifier of the device
        %   DeviceARN is the Amazon Resource Name (ARN) of the device.
        %
        %   See also quantum.backend.QuantumDeviceAWS
        DeviceARN

        %REGION - Region the device is in
        %   Region specifies the region where the device is.
        %
        %   See also quantum.backend.QuantumDeviceAWS
        Region
    end

    properties (Dependent, GetAccess = public, SetAccess = private)
        %S3Path - Path to remote storage of results
        %   S3Path is the path to a remote folder that your AWS account has
        %   access to. The result of a quantum task will be stored in this
        %   folder and then retrieved by MATLAB.
        %
        %   The S3Path must point to a folder in an S3 bucket, the name of
        %   that S3 bucket must start with "amazon-braket-" and this bucket
        %   must be accessible from your AWS_DEFAULT_REGION.
        %
        %   Default: "s3://amazon-braket-mathworks/default"
        %
        %   See also quantum.backend.QuantumDeviceAWS
        S3Path
    end

    properties(Access = private)
        S3Bucket = "amazon-braket-mathworks"
        S3BucketPrefix = "default"
        Details % getDevice struct
    end

    methods
        function s = get.S3Path(obj)
            s = "s3://" + obj.S3Bucket + "/" + obj.S3BucketPrefix;
        end
    end

    methods
        function obj = QuantumDeviceAWS(deviceName, NameValueArgs)
            arguments
                deviceName {mustBeTextScalar}
                NameValueArgs.Region {mustBeAWSRegion}
                NameValueArgs.S3Path {mustBeTextScalar}
            end

            if ismissing(deviceName)
                obj.Name = string(NaN);
                obj.DeviceARN = string(NaN);
                obj.Region = string(NaN);
                return
            end

            % Turn any char inputs to string
            deviceName = string(deviceName);
            hasRegion = isfield(NameValueArgs, 'Region');
            if hasRegion
                % has only lowercase letters, numbers, and 2 or more dashes
                NameValueArgs.Region = lower(string(NameValueArgs.Region));
            end
            userSetS3Path = isfield(NameValueArgs, 'S3Path');
            if userSetS3Path
                NameValueArgs.S3Path = string(NameValueArgs.S3Path);
            end

            % Get list of device names, ARNs and regions. If Region NVP was
            % specified, only get these specified regions.

            if hasRegion
                [devNames, devARNs, devRegions, devStatus] = ...
                    quantum.internal.aws.quickListDevices(string(NameValueArgs.Region));
            else
                % Query all regions that currently have any quantum computing devices
                [devNames, devARNs, devRegions, devStatus] = ...
                    quantum.internal.aws.quickListDevices();
            end

            % Safeguard against getting duplicate listings from aws
            % interface, just in case:
            devInfo = unique([devNames devARNs devRegions devStatus], 'rows');
            devNames = devInfo(:, 1);
            devARNs = devInfo(:, 2);
            devRegions = devInfo(:, 3);
            devStatus = devInfo(:, 4);

            % Find index of all devices that match the entered deviceName
            if startsWith(deviceName, "arn:aws:braket:")
                % Input looks like a device ARN
                ind = find(deviceName == devARNs);
            else
                ind = find(strcmpi(deviceName, devNames));
            end

            % Case no match was found
            if isempty(ind)
                if hasRegion
                    error(message("quantum:QuantumDeviceAWS:DeviceNotFoundRegion"))
                else
                    error(message("quantum:QuantumDeviceAWS:DeviceNotFound"))
                end
            end

            % Case several matches were found and region was not specified:
            % Use default region environment variable, see if this leads to
            % a unique match.
            if ~isscalar(ind) && ~hasRegion
                defaultRegion = getenv('AWS_DEFAULT_REGION');
                if ~isempty(defaultRegion) && ismember(defaultRegion, devRegions(ind))
                    [~, ii] = ismember(defaultRegion, devRegions(ind));
                    ind = ind(ii);
                end
            end

            % No unique match was found.
            if ~isscalar(ind)
                error(message("quantum:QuantumDeviceAWS:DeviceNotUnique"))
            end

            if ~strcmpi(devStatus(ind), 'online')
                error(message("quantum:QuantumDeviceAWS:DeviceNotOnline"))
            end

            obj.Name = devNames(ind);
            obj.DeviceARN = devARNs(ind);
            obj.Region = devRegions(ind);

            % Get details about device and make sure its supported
            caps = getCapabilities(obj);
            if ~isfield(caps, 'action') || ...
                    ~matches("braket_ir_openqasm_program", fields(caps.action))
                error(message("quantum:QuantumDeviceAWS:DeviceNotSupported"))
            end

            % Validate S3Path has the expected form
            if userSetS3Path
                s3p = NameValueArgs.S3Path;
                % Error if obviously invalid, then split S3Path into bucket and prefix
                if ~startsWith(s3p, "s3://amazon-braket-")
                    error(message("quantum:QuantumDeviceAWS:S3PathMustStartWith"))
                end

                slashes = strfind(s3p, '/');
                if length(slashes) < 3 || strlength(s3p) == slashes(3)
                    error(message("quantum:QuantumDeviceAWS:S3PathMustSpecifyFolder"))
                end

                % Remove one trailing "/" or "\"
                if endsWith(s3p, ("/" | "\"))
                    s3p = erase(s3p, ("/" | "\")+lineBoundary("end"));
                end

                obj.S3Bucket = extractBetween(s3p, slashes(2)+1, slashes(3)-1);
                obj.S3BucketPrefix = extractAfter(s3p, slashes(3));
            end

            % Check the full S3Path exists. If not, check bucket and create
            % folder if its found, otherwise error/warn
            if ~isfolder(obj.S3Path)
                % Need slash here to recognise bucket as a folder
                if ~isfolder("s3://"+obj.S3Bucket+"/")
                    if userSetS3Path
                        error(message('quantum:QuantumDeviceAWS:S3PathNotFound', obj.S3Path))
                    else
                        warning(message('quantum:QuantumDeviceAWS:S3PathNotFound', obj.S3Path))
                    end
                else
                    mkdir(obj.S3Path)
                end
            end
        end

        function details = fetchDetails(obj)
            %FETCHDETAILS Retrieve details about the quantum device
            %
            %   s = FETCHDETAILS(dev) returns a struct containing
            %   information about the device.
            %
            %   The output of this function may change over time. Calling
            %   FETCHDETAILS always retrieves the latest version provided
            %   by AWS.
            %
            %   See also quantum.backend.QuantumDeviceAWS.

            arguments
                obj quantum.backend.QuantumDeviceAWS
            end

            if ismissing(obj.Name)
                error(message("quantum:QuantumDevice:MissingNotSupported"))
            end

            details = getDetails(obj, true);
        end
    end

    methods (Hidden)
        function tf = isGateBased(obj)
            caps = getCapabilities(obj);
            tf = matches("braket_ir_openqasm_program", fields(caps.action));
        end

        function tf = isAnnealing(obj)
            caps = getCapabilities(obj);
            tf = matches("braket_ir_annealing_problem", fields(caps.action));
        end

        function gates = getSupportedGates(obj)
            if ismissing(obj.Name)
                error(message('quantum:QuantumDevice:MissingNotSupported'));
            end
            caps = getCapabilities(obj);
            deviceGates = caps.action.braket_ir_openqasm_program.supportedOperations;

            % This represents the supported subset of AWS Braket gates and
            % maps to the equivalent gates in MATLAB.
            AvailableGatesBraket = dictionary( ...
                ["x","y","z","rx","ry","rz","h","swap","s","t","cy","cz","ccnot","cnot", "i", "si", "ti","cphaseshift", "phaseshift", "xx", "yy", "zz"], ...
                ["x","y","z","rx","ry","rz","h","swap","s","t","cy","cz", "ccx","cx", "id", "si", "ti","cr1", "r1", "rxx", "ryy", "rzz"]);

            % Filter the device gates to those available in MATLAB and map
            % the MATLAB name to the AWS Braket name.
            available = isKey(AvailableGatesBraket, deviceGates);
            supportedBraket = string(deviceGates(available));
            gates = dictionary(AvailableGatesBraket(supportedBraket), supportedBraket);
        end

        function qasm = generateQASM(obj, circuit)
            % Generate OpenQASM for the AWS device
            gates = obj.getSupportedGates();
            qasm = generateQASM(circuit, Unpack=true, SupportedGates=gates);
            qasm = erase(qasm, 'include "stdgates.inc";');
        end

        function task = sendQASM(obj, qasm, NameValueArgs)
            if ismissing(obj.Name)
                error(message('quantum:QuantumDevice:MissingNotSupported'));
            end

            action.braketSchemaHeader = struct(name="braket.ir.openqasm.program", version="1");

            action.source = qasm;
            action = jsonencode(action);

            task = send(obj, action, NameValueArgs.NumShots);
        end

        function task = sendCircuit(obj, circuit, NameValueArgs)
            arguments
                obj quantum.backend.QuantumDeviceAWS
                circuit quantumCircuit
                NameValueArgs.NumShots (1,1) {mustBeInteger, mustBePositive} = 100
                NameValueArgs.Observable
            end

            if ~isscalar(circuit)
                % AWS only supports one circuit per task
                error(message('quantum:QuantumDeviceAWS:invalidCircuit'))
            end

            qasm = generateQASM(obj, circuit);

            if isfield(NameValueArgs, "Observable")
                obs = NameValueArgs.Observable;
                if ~(isa(obs, "observable") && isscalar(obs) && isscalar(obs.Paulis))
                    % AWS only supports one observable per task
                    error(message('quantum:QuantumDeviceAWS:invalidObservable'));
                end
                if circuit.NumQubits~=obs.NumQubits
                    error(message("quantum:observable:observeIncorrectNumQubits"))
                end
                N = obs.NumQubits;
                % Result types cannot be used with classical registers
                qasm = erase(qasm, "bit["+N+"] c;");
                % Remove measure because this implies samples of Z observable
                qasm = erase(qasm, 'c = measure q;');
                % Add observable weight as comment to be multiplied by
                % task. Weights are not directly supported.
                w = obs.Weights;
                if isa(w, 'single')
                    wStr = num2str(w, 9);
                else
                    wStr = num2str(w, 17);
                end
                wComment = newline+"//Start weight"+newline+"//"+wStr+newline+"//End weight";
                qasm = insertAfter(qasm, 'OPENQASM 3.0;', wComment);
                % Request expectation result type from the server
                p = lower(char(obs.Paulis)') + "(q["+string(0:N-1).'+"])";
                qasm = qasm + newline + "#pragma braket result expectation "+join(p, " @ ");
            end

            task = sendQASM(obj, qasm, NameValueArgs);
        end

        function s = saveobj(~)
            s = struct;
            warning(message("quantum:QuantumDevice:SaveNotSupported"));
        end
    end

    methods (Hidden, Static)
        function obj = loadobj(~)
            warning(message('quantum:QuantumDevice:LoadNotSupported'));
            obj = quantum.backend.QuantumDeviceAWS(string(NaN));
        end
    end

    methods (Access = private)
        function details = getDetails(obj, forceQuery)
            if ismissing(obj.Name)
                error(message('quantum:QuantumDevice:MissingNotSupported'));
            end

            if isempty(obj.Details) || forceQuery
                obj.Details = quantum.internal.aws.getDevice(obj.DeviceARN, obj.Region);
            end
            details = obj.Details;
        end

        function capabilities = getCapabilities(obj)
            %NOTE: deviceCapabilities field is not expected to change often
            details = getDetails(obj, false);
            capabilities = jsondecode(details.deviceCapabilities);
        end

        function task = send(device, action, numShots)
            % Use for both gate-based and annealing to send a quantum task
            taskARN = quantum.internal.aws.createQuantumTask(...
                device.Region, device.DeviceARN, numShots, device.S3Bucket, ...
                device.S3BucketPrefix, action);
            task = quantum.backend.QuantumTaskAWS(taskARN);
        end
    end
end

function mustBeAWSRegion(str)
mustBeTextScalar(str)
dashes = strfind(str,"-");
if ~matches(str, regexpPattern('[a-zA-Z0-9-]+')) || length(dashes)<2
    error(message('quantum:QuantumDeviceAWS:InvalidRegion'));
end
end
