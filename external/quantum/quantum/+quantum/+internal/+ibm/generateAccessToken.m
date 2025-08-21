function credentials = generateAccessToken(credentials)
% Internal use only.

% Copyright 2022-2024 The MathWorks, Inc.

if isequal(credentials.channel, "ibm_cloud")

    % The Cloud channel uses the permanent token to generated a temporary
    % access token used for authentication. This is generated whenever a device or
    % task is directly constructed.
    url = "https://iam.cloud.ibm.com/identity/token";
    authHeader = matlab.net.http.HeaderField("Content-Type","application/x-www-form-urlencoded");
    body = "grant_type=urn:ibm:params:oauth:grant-type:apikey&apikey="+credentials.token;

    method = matlab.net.http.RequestMethod.POST;
    request = matlab.net.http.RequestMessage(method, authHeader, body);

    response = quantum.internal.ibm.sendRequest(request, url);
    credentials.accessToken = string(response.access_token);
end
