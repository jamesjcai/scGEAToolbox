function [done, res] = callOpenAIChatNVIDIA(prompt, model)

done = false;
res = [];

if isempty(which('openAIChat'))
    error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
end

if nargin < 2, model = "minimaxai/minimax-m2.5"; end
if nargin < 1, prompt = 'Why is the sky blue?'; end


preftagname = 'llapikeyenvfile';
apikeyfile = getpref('scgeatoolbox', preftagname);
loadenv(apikeyfile,"FileType","env");
apikey = getenv("NVIDIA_API_KEY");
apibase = getenv("NVIDIA_API_BASE");

if isempty(apibase)
    apibase="https://integrate.api.nvidia.com/v1";
end

fprintf('Sending request to OpenAI Chat API...\n');
try
    chat = openAIChat("", APIKey=apikey, ...
        ModelName=model, TimeOut=1200, BaseURL=apibase);
    res = chat.generate(prompt);
    done = true;
catch ME
    fprintf('Error in chat completion: %s\n', ME.message);
    return;
end
if done
    fprintf('Response received successfully.\n');
end
% disp(res);
end
