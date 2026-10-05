function values = i_readkeyfile(keyfile)
%I_READKEYFILE The NAME=VALUE pairs of an llm_api_key.env file.
%   values = llm.i_readkeyfile() reads the file named by the scgeatoolbox
%   'llapikeyenvfile' preference.
%   values = llm.i_readkeyfile(keyfile) reads KEYFILE instead.
%
%   VALUES is a string-to-string dictionary, empty when there is no file.
%   Nothing is copied into the environment. LOADENV with no output sets
%   every pair as an environment variable, so a key in the file shadowed
%   the same key stored in the MATLAB vault; see llm.i_getapikey.
%
%   See also llm.i_getapikey, loadenv.

arguments
    keyfile {mustBeTextScalar} = ""
end

values = dictionary(string.empty, string.empty);
keyfile = string(keyfile);
if strlength(keyfile) == 0 && ispref("scgeatoolbox", "llapikeyenvfile")
    keyfile = string(getpref("scgeatoolbox", "llapikeyenvfile"));
end
if strlength(keyfile) == 0 || ~isfile(keyfile)
    return;
end
values = loadenv(keyfile, FileType="env");
end
