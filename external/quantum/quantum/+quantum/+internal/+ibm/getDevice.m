function details = getDevice(credentials, deviceName)
% Internal use only.

% There are 4 different endpoints of information available that each get
% their own field.

% Copyright 2022-2024 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);
rootURL = rootURL+"/backends/"+deviceName;

method = matlab.net.http.RequestMethod.GET;
request = matlab.net.http.RequestMessage(method, authHeader);

% Supported for all devices
details.Status = quantum.internal.ibm.sendRequest(request, rootURL+"/status");
details.Configuration = quantum.internal.ibm.sendRequest(request, rootURL+"/configuration");

% Supported for real hardware
props = quantum.internal.ibm.sendRequest(request, rootURL+"/properties");
defaults = quantum.internal.ibm.sendRequest(request, rootURL+"/defaults");

% Ensure properties and defaults are non-empty structs to safeguard against
% users who may have access to simulators, which would return an empty array.
if isstruct(props) && ~isempty(fieldnames(props))
    details.Properties = props;
end
if isstruct(defaults) && ~isempty(fieldnames(defaults))
    details.Defaults = defaults;
end
end