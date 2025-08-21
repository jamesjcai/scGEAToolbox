function [Name, DeviceARN, Region, Status, Provider, deviceType] = quickListDevices(regions)
% Helper for QuantumDeviceAWS

%   Copyright 2022-2024 The MathWorks, Inc.

if nargin == 0
    % All regions where devices are available 
    regions = ["us-east-1", "us-west-1", "us-west-2", "eu-west-2", "eu-north-1"];
end

regions = string(regions);

Name = string.empty;
Status = string.empty;
DeviceARN = string.empty;
Provider = string.empty;
Region = string.empty;
deviceType = string.empty;

numDev = 0;
for r = 1:length(regions)
    devs = quantum.internal.aws.searchDevices(regions(r));
    % add region field to override the default in get-devices
    for d = 1:length(devs)
        if ~strcmpi(devs(d).deviceStatus, 'retired')
            numDev = numDev+1;
            Name(numDev,1) = devs(d).deviceName;
            Status(numDev,1) = devs(d).deviceStatus;
            DeviceARN(numDev,1) = devs(d).deviceArn;
            Provider(numDev,1) = devs(d).providerName;
            deviceType(numDev,1) = devs(d).deviceType;
            Region(numDev,1) = regions(r);
        end
    end
end
