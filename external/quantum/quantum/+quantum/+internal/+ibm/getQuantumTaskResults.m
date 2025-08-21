function results = getQuantumTaskResults(credentials, taskID)
% Internal use only. Returns struct of results associated with the taskID.

% Copyright 2022-2024 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);

method = matlab.net.http.RequestMethod.GET;
request = matlab.net.http.RequestMessage(method, authHeader);

response = quantum.internal.ibm.sendRequest(request, rootURL+"/jobs/"+taskID+"/results");

% Safeguard against server returning struct instead of JSON string.
if isstring(response)
    results = jsondecode(response);
else
    assert(isstruct(response))
    results = response;
end