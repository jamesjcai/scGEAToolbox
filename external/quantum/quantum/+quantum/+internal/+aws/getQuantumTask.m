function taskDetails = getQuantumTask(taskARN, taskRegion)
%GETQUANTUMTASK Query task details

%   Copyright 2022-2023 The MathWorks, Inc.

userAgentString = quantum.internal.aws.userAgentString();
quantum.internal.aws.addDLLLocationToSystemPathIfNecessary();
taskDetails = quantum.internal.aws.getQuantumTaskMex(taskARN, taskRegion, userAgentString);

