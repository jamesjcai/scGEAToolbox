function taskARN = createQuantumTask(region, deviceARN, numShots, s3Bucket, s3BucketPrefix, action)
%CREATEQUANTUMTASK Send a new quantum task to aws

%   Copyright 2022-2023 The MathWorks, Inc.

userAgentString = quantum.internal.aws.userAgentString();
quantum.internal.aws.addDLLLocationToSystemPathIfNecessary();
taskARN = quantum.internal.aws.createQuantumTaskMex(...
  region, deviceARN, numShots, s3Bucket, s3BucketPrefix, action, userAgentString);

