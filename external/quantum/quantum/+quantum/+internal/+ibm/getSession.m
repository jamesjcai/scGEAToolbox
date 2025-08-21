function details = getSession(credentials, sessionID)
% Internal use only.

% Copyright 2022-2024 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);

method = matlab.net.http.RequestMethod.GET;
request = matlab.net.http.RequestMessage(method, authHeader);

details = quantum.internal.ibm.sendRequest(request, rootURL+"/sessions/"+sessionID);
end