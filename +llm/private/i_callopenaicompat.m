function [done, res] = i_callopenaicompat(prompt, model, keyVars, baseVar, defaultBase, label)
%I_CALLOPENAICOMPAT One prompt to an OpenAI-compatible chat endpoint.
%   [done, res] = i_callopenaicompat(prompt, model, keyVars, baseVar, ...
%       defaultBase, label)
%
%   The shared body of llm.callOpenAIChat, llm.callOpenAIChatTAMU and
%   llm.callOpenAIChatNVIDIA, which were one function written three times
%   apart from the variable names (the TAMU copy also printed the reply
%   twice). KEYVARS are environment variables tried in order for the API
%   key; BASEVAR names the base-URL variable, DEFAULTBASE is used when it is
%   unset ("" for the OpenAI default endpoint). The env file named by the
%   scgeatoolbox 'llapikeyenvfile' preference is loaded first when there is
%   one; the copies called GETPREF unconditionally and threw when it was
%   unset, instead of trying the environment as it stood.
%
%   done is true on a reply; on failure res is [] and the reason is printed
%   as "Error in chat completion: ..." (llm.i_askllm passes that on).

done = false;
res = [];

if isempty(which('openAIChat'))
    error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
end

if ispref('scgeatoolbox', 'llapikeyenvfile')
    apikeyfile = getpref('scgeatoolbox', 'llapikeyenvfile');
    if ~isempty(apikeyfile) && isfile(apikeyfile)
        loadenv(apikeyfile, "FileType", "env");
    end
end

apikey = "";
for v = string(keyVars)
    apikey = string(getenv(v));
    if strlength(apikey) > 0, break; end
end
if strlength(apikey) == 0
    fprintf('Error in chat completion: %s is not set.\n', strjoin(string(keyVars), " or "));
    return;
end

apibase = "";
if strlength(string(baseVar)) > 0
    apibase = string(getenv(baseVar));
end
if strlength(apibase) == 0
    apibase = string(defaultBase);
end

fprintf('Sending request to %s Chat API...\n', label);
try
    if strlength(apibase) > 0
        chat = openAIChat("", APIKey=apikey, ModelName=model, TimeOut=1200, BaseURL=apibase);
    else
        chat = openAIChat("", APIKey=apikey, ModelName=model, TimeOut=1200);
    end
    res = chat.generate(prompt);
    done = true;
catch ME
    fprintf('Error in chat completion: %s\n', ME.message);
    return;
end
fprintf('Response received successfully.\n');
end
