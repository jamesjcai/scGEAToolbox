function [done, tested] = i_checkllm(apikeyfile, provider, ~)
% The third argument, the parent figure, is unused since the env-file
% picker moved to gui.i_setllmmodel; it stays for existing callers.

done = false;

if nargin<2, provider = []; end
if nargin<1, apikeyfile = []; end

% No env file is required: keys can also come from the environment or the
% MATLAB vault (llm.i_getapikey). With no file and no preference, those two
% are all that is searched.
if isempty(apikeyfile) && ispref('scgeatoolbox', 'llapikeyenvfile')
    apikeyfile = getpref('scgeatoolbox', 'llapikeyenvfile');
end
if isempty(apikeyfile), apikeyfile = ""; end
fileValues = llm.i_readkeyfile(apikeyfile);
if strlength(apikeyfile) > 0
    fprintf('Using apikeyfile %s\n', apikeyfile);
end


preftagname = 'llmodelprovider';
s = getpref('scgeatoolbox', preftagname);
providermodel = strsplit(s,':');

if isempty(provider)
        provider = providermodel{1};
    end

fprintf('Using LLM provider: %s\n', provider);

model = strjoin(providermodel(2:end), ':');
fprintf('Using LLM model: %s\n', model);

prompt = "What model are you?";

% The switch below sends a test prompt for these providers only. Any other
% falls through to DONE = TRUE untested, which the caller used to report
% as a success; TESTED lets it say so instead.
tested = ismember(string(provider), ["Ollama", "OpenAI", "TAMUAIChat", "Gemini"]);

switch provider
        case 'Ollama'
            try
                chat = ollamaChat(model, TimeOut = 1200);
                feedbk = chat.generate(prompt);
            catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                return;
            end
            disp(feedbk);
        case 'OpenAI'
            apikey = llm.i_getapikey(["OPENAI_API_KEY", "OpenAI_API_KEY"], fileValues);
            if isempty(apikey), return; end
            try
                chat = openAIChat("", APIKey=apikey, ...
                    ModelName=model, TimeOut=1200);
                feedbk = chat.generate(prompt);
            catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                return;
            end
            disp(feedbk);
        case 'TAMUAIChat'
            OPEN_WEBUI_API_ENDPOINT = "https://chat-api.tamu.ai";
            OPEN_WEBUI_API_KEY = llm.i_getapikey("TAMUAI_API_KEY", fileValues);

            if isempty(OPEN_WEBUI_API_KEY), return; end
            chat_url = sprintf('%s/api/chat/completions', OPEN_WEBUI_API_ENDPOINT);

            % Create request body structure
            messages_cell = {struct('role', 'user', 'content', prompt)};
            body_struct = struct('model', model, ...
                                'stream', false, ...
                                'messages', {messages_cell});
            % Debug: Display the request body structure
            fprintf('Request body structure:\n');
            disp(body_struct);

            % Set up options for webwrite (POST request)
            post_options = weboptions('MediaType', 'application/json', ...
                                     'RequestMethod', 'POST', ...
                                     'HeaderFields', {'Authorization', ...
                                     sprintf('Bearer %s', OPEN_WEBUI_API_KEY)});

            % Make the chat completion request
            try
                chat_response = webwrite(chat_url, body_struct, post_options);
                % Display chat response as JSON
                chat_json = jsonencode(chat_response);
                fprintf('Chat completion response:\n%s\n', chat_json);

            catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                return;
            end
        case 'Gemini'
            % CALLGEMINI takes the key FILE and reads the key from it. This
            % passed the key itself, which then failed ISFILE, so the key
            % was never loaded and the Gemini check always failed.
            try
                [ok, response] = llm.callGemini(apikeyfile, prompt, model);
                if ~ok
                    fprintf('Gemini returned an error.\n');
                    return;
                end
             catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                return;
            end
            disp(response);
    end
done = true;
end
