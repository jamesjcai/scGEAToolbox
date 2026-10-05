function [key, source] = i_getapikey(names, fileValues)
%I_GETAPIKEY The first API key found under any of NAMES.
%   key = llm.i_getapikey(names)
%   key = llm.i_getapikey(names, fileValues)
%   [key, source] = llm.i_getapikey(___)
%
%   NAMES is one variable name or several tried in order, for example
%   ["OPENAI_API_KEY", "OpenAI_API_KEY"]. Each place is searched for every
%   name before the next place is tried:
%
%     1. the environment (GETENV)
%     2. the MATLAB vault (GETSECRET), where llm.i_storeapikey puts keys
%     3. FILEVALUES, the pairs read by llm.i_readkeyfile; by default the
%        file named by the 'llapikeyenvfile' preference. Pass an empty
%        dictionary to skip the file.
%
%   The environment comes first, as it does for secrets in a deployed
%   application, so a key set for one session or by a test still wins.
%
%   KEY is a char vector, '' when nothing was found, so ISEMPTY tests it as
%   it tested GETENV's result. SOURCE is "environment", "vault", "file" or
%   "".
%
%   Setting the environment variable SCGEA_SECRET_VAULT to "off" skips the
%   vault. The tests do that, so keys in the developer's vault do not change
%   their results.
%
%   See also llm.i_storeapikey, llm.i_readkeyfile, getSecret.

arguments
    names {mustBeText}
    fileValues = llm.i_readkeyfile()
end

names = string(names);
key = '';
source = "";

for name = names
    v = getenv(name);
    if ~isempty(v)
        key = v;
        source = "environment";
        return;
    end
end

if ~strcmpi(getenv("SCGEA_SECRET_VAULT"), "off")
    for name = names
        v = i_vaultvalue(name);
        if strlength(v) > 0
            key = char(v);
            source = "vault";
            return;
        end
    end
end

for name = names
    if isConfigured(fileValues) && isKey(fileValues, name)
        v = fileValues(name);
        if strlength(v) > 0
            key = char(v);
            source = "file";
            return;
        end
    end
end
end

function v = i_vaultvalue(name)
% The vault value of NAME, or "" when there is none.
v = "";
try
    if isSecret(name)
        v = string(getSecret(name));
    end
catch ME
    % A vault that cannot be opened (a locked keychain, a headless
    % session) leaves the env file to supply the key.
    fprintf("MATLAB vault unavailable (%s); trying the env file.\n", ME.message);
end
end
