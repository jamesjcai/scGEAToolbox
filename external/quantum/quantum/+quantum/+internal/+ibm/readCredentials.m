function [credentials, accountname] = readCredentials(fileName, accountName)
% Internal use only. Returns struct of credentials from the JSON file. This
% does not authenticate with the IBM Runtime service.

% Copyright 2022-2025 The MathWorks, Inc.

% Read content of JSON file into struct
content = fileread(fileName);
accountsStruct = jsondecode(content);

if ~(isscalar(accountsStruct) && isstruct(accountsStruct) && ~isempty(fieldnames(accountsStruct)))
    error(message("quantum:QuantumDeviceIBM:InvalidFileContent"))
end

accountNames = fieldnames(accountsStruct);

% Get struct of account credentials
if strcmp(accountName, "")
    % Use first account by default
    accountname = string(accountNames{1});
    account = accountsStruct.(accountname);
else
    accountname = matlab.lang.makeValidName(accountName);
    if ~isfield(accountsStruct, accountname)
        error(message("quantum:QuantumDeviceIBM:AccountNotFound"))
    end
    account = accountsStruct.(accountname);
end

if ~(isscalar(account) && isstruct(account) && length(fieldnames(account)) >= 3)
    error(message("quantum:QuantumDeviceIBM:InvalidCredentialsContent", accountname))
end

names = fieldnames(account);
values = struct2cell(account);

% Locate case-insensitive matches for credentials
creds = ["channel", "instance", "token"];
credentials = struct;
for c = creds
    tf = matches(names,c, IgnoreCase=true);
    if nnz(tf)~=1
        error(message("quantum:QuantumDeviceIBM:InvalidCredentialsContent", accountname))
    end
    credentials.(c) = string(values{tf});
end

if credentials.channel~="ibm_cloud"
    error(message("quantum:QuantumDeviceIBM:InvalidChannel"))
elseif ~contains(credentials.instance, "crn:v1:bluemix:public:quantum-computing")
    error(message("quantum:QuantumDeviceIBM:InvalidCloudInstance"))
end

quantum.internal.ibm.showLegalDisclaimer()

end