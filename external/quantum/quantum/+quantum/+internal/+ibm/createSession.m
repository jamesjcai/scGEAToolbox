function sessionID = createSession(credentials, sessionInfo)
% Internal use only.

% Copyright 2022-2025 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);

method = matlab.net.http.RequestMethod.POST;
request = matlab.net.http.RequestMessage(method, authHeader, sessionInfo);

response = quantum.internal.ibm.sendRequest(request, rootURL+"/sessions");
sessionID = string(response.id);
end