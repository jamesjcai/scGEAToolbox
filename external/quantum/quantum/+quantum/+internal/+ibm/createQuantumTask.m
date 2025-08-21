function taskID = createQuantumTask(credentials, taskInfo)
% Internal use only.

% Copyright 2022-2024 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);

method = matlab.net.http.RequestMethod.POST;
request = matlab.net.http.RequestMessage(method, authHeader, taskInfo);

response = quantum.internal.ibm.sendRequest(request, rootURL+"/jobs");
taskID = string(response.id);
end