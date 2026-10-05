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
%   unset ("" for the OpenAI default endpoint). The key is looked up by
%   llm.i_getapikey (environment, MATLAB vault, then the env file named by
%   the 'llapikeyenvfile' preference); the base URL in the environment, then
%   the env file.
%
%   done is true on a reply; on failure res is [] and the reason is printed
%   as "Error in chat completion: ..." (llm.i_askllm passes that on).

done = false;
res = [];

if isempty(which('openAIChat'))
    error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
end

fileValues = llm.i_readkeyfile();
apikey = string(llm.i_getapikey(keyVars, fileValues));
if strlength(apikey) == 0
    fprintf(['Error in chat completion: %s is not set. Store it with ' ...
        'llm.i_storeapikey("%s") or add it to llm_api_key.env.\n'], ...
        strjoin(string(keyVars), " or "), string(keyVars(1)));
    return;
end

apibase = "";
baseVar = string(baseVar);
if strlength(baseVar) > 0
    apibase = string(getenv(baseVar));
    if strlength(apibase) == 0 && isConfigured(fileValues) && isKey(fileValues, baseVar)
        apibase = fileValues(baseVar);
    end
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
