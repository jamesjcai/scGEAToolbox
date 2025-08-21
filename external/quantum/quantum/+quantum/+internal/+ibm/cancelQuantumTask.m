function cancelQuantumTask(credentials, taskID)
% Internal use only.

% Copyright 2022-2024 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);

% The Runtime API uses a POST request to cancel but requires no body
% content, so include empty struct to avoid warning from matlab.net.http
method = matlab.net.http.RequestMethod.POST;
request = matlab.net.http.RequestMessage(method, authHeader, struct);

quantum.internal.ibm.sendRequest(request, rootURL+"/jobs/"+taskID+"/cancel");
end
