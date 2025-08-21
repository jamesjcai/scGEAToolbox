function r = devicesForTabCompletion()
% Helper for tab completion in QuantumDeviceAWS constructor

%   Copyright 2023 The MathWorks, Inc.

persistent devices timestamp;

if isempty(timestamp) || (datetime('now') - timestamp > hours(2))
    [devices, deviceARN, region] = quantum.internal.aws.quickListDevices();
    timestamp = datetime('now');
    
    % Uniquify the list by simple device names:
    [devices, ind] = unique(devices);
    deviceARN = deviceARN(ind);
    region = region(ind);

    % Check if device is supported:
    isSupported = false(size(devices));
    for ii=1:length(devices)
        if strcmpi(devices(ii), 'SV1') || strcmpi(devices(ii), 'dm1') || strcmpi(devices(ii), 'TN1')
            % Hard-code the simulators for speed, no need to look up
            % their capabilities
            isSupported(ii) = true;
            continue;
        end
           
        details = quantum.internal.aws.getDevice(deviceARN(ii), region(ii));
        capabilities = jsondecode(details.deviceCapabilities);
        isSupported(ii) = isfield(capabilities, 'action') && ...
            matches("braket_ir_openqasm_program", fields(capabilities.action));
    end

    % Only suggest supported devices
    devices = devices(isSupported);
end

r = devices;