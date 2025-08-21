function devices = listDevices(credentials)
% Internal use only. Returns string array of device names or an empty
% string array when none are available.

% Copyright 2022-2024 The MathWorks, Inc.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);

method = matlab.net.http.RequestMethod.GET;
request = matlab.net.http.RequestMessage(method, authHeader);

url = rootURL+"/backends";

if isequal(credentials.channel, "ibm_quantum")
    % Include the provider query parameter set as the instance to filter
    % devices
    url = url+"?provider="+credentials.instance;
end

response = quantum.internal.ibm.sendRequest(request, url);
% Convert cellstr to string array
devices = string(response.devices);
end