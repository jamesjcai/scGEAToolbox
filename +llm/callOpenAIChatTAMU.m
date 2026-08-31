function [done, res] = callOpenAIChatTAMU(prompt, model)

done = false;
res = [];

if isempty(which('openAIChat'))
    error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
end

if nargin < 2, model = "protected.gemini-2.0-flash-lite"; end
if nargin < 1, prompt = 'Why is the sky blue?'; end


preftagname = 'llapikeyenvfile';
apikeyfile = getpref('scgeatoolbox', preftagname);
loadenv(apikeyfile,"FileType","env");
apikey = getenv("TAMUAI_API_KEY");
apibase = getenv("TAMUAI_API_BASE");

if isempty(apibase)
    apibase = "https://chat-api.tamu.ai/api";
end

fprintf('Sending request to OpenAI Chat API...\n');
try
    %chat = openAIChat("APIKey", apikey, ...
    %    "ModelName", model,TimeOut = 1200, BaseURL = "https://chat-api.tamu.ai/api/chat/completions");
    chat = openAIChat("", APIKey=apikey, ...
        ModelName=model, TimeOut=1200, BaseURL=apibase);
    res = chat.generate(prompt);
    done = true;
catch ME
    fprintf('Error in chat completion: %s\n', ME.message);
    return;
end
if done
    disp(res);
    fprintf('Response received successfully.\n');
end
    disp(res);
end
