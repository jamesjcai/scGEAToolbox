function [key, done] = i_storeapikey(name)
%I_STOREAPIKEY Ask for an API key and keep it in the MATLAB vault.
%   key = llm.i_storeapikey(name)
%   [key, done] = llm.i_storeapikey(name)
%
%   Opens MATLAB's masked SETSECRET dialog for the secret NAME (for example
%   "TAMUAI_API_KEY"), replacing any value already stored under it. The key
%   is kept encrypted in the user's vault, so it is never written to a file
%   or a preference; llm.i_getapikey finds it in every later session.
%
%   SETSECRET takes no value argument, by design: a key typed into code
%   would land in the command history. So keys in an env file cannot be
%   moved into the vault by a script, only entered here.
%
%   KEY is the stored value as a char vector, '' when the dialog was
%   cancelled or the vault is unavailable; DONE is true when a key was
%   stored.
%
%   See also llm.i_getapikey, setSecret, removeSecret.

arguments
    name {mustBeTextScalar}
end

key = '';
done = false;
try
    setSecret(name, Overwrite=true);
catch ME
    fprintf("Could not store %s in the MATLAB vault: %s\n", name, ME.message);
    return;
end
if isSecret(name)
    key = char(getSecret(name));
    done = ~isempty(key);
end
end
