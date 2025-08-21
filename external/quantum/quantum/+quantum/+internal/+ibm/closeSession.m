function closeSession(credentials, sessionID)
% Internal use only.

% This is just used by tests to verify behavior. Sessions are automatically
% closed by the server.

% Copyright 2022-2024 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);

method = matlab.net.http.RequestMethod.DELETE;
request = matlab.net.http.RequestMessage(method, authHeader);

quantum.internal.ibm.sendRequest(request, rootURL+"/sessions/"+sessionID+"/close");
end