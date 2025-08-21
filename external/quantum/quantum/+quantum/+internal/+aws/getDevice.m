function details = getDevice(deviceARN, region)
%GETDEVICE Query device details

%   Copyright 2022-2023 The MathWorks, Inc.

userAgentString = quantum.internal.aws.userAgentString();
quantum.internal.aws.addDLLLocationToSystemPathIfNecessary();
details = quantum.internal.aws.getDeviceMex(deviceARN, region, userAgentString);
