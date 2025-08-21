function [rootURL, authHeader] = connectionInfo(credentials)
% Internal use only. Static data used for connecting to the Runtime Service

% Copyright 2022-2024 The MathWorks, Inc.

% The Qiskit version string was provided by IBM to identify non-Python
% users. This does not need to be updated.
userAgent = "qiskit-version-2/0.39.2/"+quantum.internal.versionString();

if isequal(credentials.channel, "ibm_quantum")
    authHeader = [
        matlab.net.http.HeaderField("Authorization", "Bearer "+credentials.token), ...
        matlab.net.http.HeaderField("x-qx-client-application", userAgent)
        ];

    rootURL = "https://api.quantum-computing.ibm.com/runtime";
else
    authHeader = [
        matlab.net.http.HeaderField("Authorization", "Bearer "+credentials.accessToken), ...
        matlab.net.http.HeaderField("Service-CRN", credentials.instance), ...
        matlab.net.http.HeaderField("x-qx-client-application", userAgent)
        ];

    rootURL = "https://us-east.quantum-computing.cloud.ibm.com";
end
end

