function [done, res] = callOllama(prompt, model)

done = false;

    if isempty(which('ollamaChat'))
        error('Needs the Add-On of Large Language Models (LLMs) with MATLAB');
    end

    % https://www.mathworks.com/matlabcentral/fileexchange/163796-large-language-models-llms-with-matlab/
    % Add-On "Large Language Models (LLMs) with MATLAB".

    if nargin < 2, model = "llama3.2"; end
    if nargin < 1, prompt = 'Why is the sky blue?'; end

    fprintf('Sending request to Ollama Chat API...\n');
    try
        chat = ollamaChat(model, TimeOut = 1200);
        res = generate(chat, prompt);
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