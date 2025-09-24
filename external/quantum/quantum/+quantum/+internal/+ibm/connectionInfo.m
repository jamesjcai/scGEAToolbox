function [rootURL, authHeader] = connectionInfo(credentials)
% Internal use only. Static data used for connecting to the Runtime Service

% Copyright 2022-2025 The MathWorks, Inc.

% The Qiskit version string was provided by IBM to identify non-Python
% users. This does not need to be updated.
userAgent = "qiskit-version-2/0.39.2/"+quantum.internal.versionString();

authHeader = [
    matlab.net.http.HeaderField("Authorization", "Bearer "+credentials.accessToken), ...
    matlab.net.http.HeaderField("Service-CRN", credentials.instance), ...
    matlab.net.http.HeaderField("x-qx-client-application", userAgent)
    ];

% These are the only two regions that exist and covers all users.
inst = credentials.instance;
if contains(inst, "us-east", IgnoreCase=true)
    rootURL = "https://quantum.cloud.ibm.com/api/v1";
elseif contains(inst, "eu-de", IgnoreCase=true)
    rootURL = "https://eu-de.quantum.cloud.ibm.com/api/v1";
else
    error(message("quantum:QuantumDeviceIBM:InvalidCloudInstance"))
end

end

