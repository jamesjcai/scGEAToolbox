function cals = getDeviceCalibration(credentials, deviceName)
% Internal use only. Returns 2x2xN array of measurement error for N qubits.
% Each 2x2 matrix represents probability of correctly measuring 0 and 1.

% Copyright 2024 The MathWorks, Inc.

% References:
% [1] P. D. Nation, H. Kang, N. Sundaresan, J. M. Gambetta. "Scalable
% Mitigation of Measurement Errors on Quantum Computers." PRX Quantum,
% Nov. 2021. American Physical Society.

[rootURL, authHeader] = quantum.internal.ibm.connectionInfo(credentials);
rootURL = rootURL+"/backends/"+deviceName;

method = matlab.net.http.RequestMethod.GET;
request = matlab.net.http.RequestMessage(method, authHeader);

props = quantum.internal.ibm.sendRequest(request, rootURL+"/properties");

% N-by-8 struct with various metrics.
info = props.qubits;
N = size(info,1);
cals = zeros(2,2, N);
for k = 1:N
    tbl = struct2table(info(k,:));
    idx01 = matches(tbl.name, "prob_meas0_prep1");
    idx10 = matches(tbl.name, "prob_meas1_prep0");
    p01 = tbl(idx01,:).value;
    p10 = tbl(idx10,:).value;
    % Equation 3 [1]
    cals(:,:,k) = [1-p10 p01; p10 1-p01];
end
end