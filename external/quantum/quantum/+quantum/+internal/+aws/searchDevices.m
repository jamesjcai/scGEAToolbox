function devices = searchDevices(region)
%SEARCHDEVICES Query all devices available in region

%   Copyright 2022-2023 The MathWorks, Inc.

userAgentString = quantum.internal.aws.userAgentString();
quantum.internal.aws.addDLLLocationToSystemPathIfNecessary();
devices = quantum.internal.aws.searchDevicesMex(region, userAgentString);

