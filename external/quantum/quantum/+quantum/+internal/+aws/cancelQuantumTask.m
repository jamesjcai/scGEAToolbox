function cancelQuantumTask(taskARN, taskRegion)
%CANCELQUANTUMTASK Cancel an existing quantum task

%   Copyright 2022-2023 The MathWorks, Inc.

userAgentString = quantum.internal.aws.userAgentString();
quantum.internal.aws.addDLLLocationToSystemPathIfNecessary();
status = quantum.internal.aws.cancelQuantumTaskMex(taskARN, taskRegion, userAgentString);

