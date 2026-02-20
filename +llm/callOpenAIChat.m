function [done, res] = callOpenAIChat(prompt, model)

    done = false;
    res = [];

    if isempty(which('openAIChat'))
        error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
    end

    if nargin < 2, model = "gpt-4.1"; end
    if nargin < 1, prompt = 'Why is the sky blue?'; end


    preftagname = 'llapikeyenvfile';
    apikeyfile = getpref('scgeatoolbox', preftagname);
    loadenv(apikeyfile,"FileType","env");
    apikey = getenv("OPENAI_API_KEY");
    
    fprintf('Sending request to OpenAI Chat API...\n');
    try
        chat = openAIChat("APIKey", apikey, ...
            "ModelName", model,TimeOut = 1200);
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



