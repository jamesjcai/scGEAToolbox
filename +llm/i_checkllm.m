function [done] = i_checkllm(apikeyfile, provider, parentfig)

done = false;

if nargin<3, parentfig = []; end
if nargin<2, provider = []; end
if nargin<1, apikeyfile = []; end

if isempty(apikeyfile)
    preftagname = 'llapikeyenvfile';
    if ~ispref('scgeatoolbox', preftagname)
        if ~strcmp('Yes', gui.myQuestdlg(parentfig, 'Locate llm_api_key.env?')), return; end
        [file, path] = uigetfile('llm_api_key.env', 'Select File');
        if isequal(file, 0), return; end
        apikeyfile = fullfile(path, file);
        setpref('scgeatoolbox', preftagname, apikeyfile);
        gui.myHelpdlg(parentfig, "llm_api_key.env is located successfully.");
    else
        apikeyfile = getpref('scgeatoolbox', preftagname);
    end
end
fprintf('Using apikeyfile %s\n', apikeyfile);


    preftagname = 'llmodelprovider';
    s = getpref('scgeatoolbox', preftagname);
    providermodel = strsplit(s,':');

    if isempty(provider)
        provider = providermodel{1};
    end

fprintf('Using LLM provider: %s\n', provider);

model = providermodel{2};
fprintf('Using LLM model: %s\n', model);

prompt = "What model are you?";
    
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
            loadenv(apikeyfile, "FileType", "env");
            apikey = getenv("OPENAI_API_KEY");
            if isempty(apikey), return; end
            try
                chat = openAIChat("APIKey",apikey, ...
                    "ModelName", model, TimeOut = 1200);                
                feedbk = chat.generate(prompt);
            catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                return;
            end
            disp(feedbk);        
        case 'TAMUAIChat'

            loadenv(apikeyfile, "FileType", "env");
            OPEN_WEBUI_API_ENDPOINT = "https://chat-api.tamu.ai";
            OPEN_WEBUI_API_KEY = getenv("TAMUAI_API_KEY");

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
            loadenv(apikeyfile,"FileType","env");
            apiKey = getenv("GEMINI_API_KEY");            
            try
                response = llm.callGemini(apiKey, prompt, model);
                % response = llm.geminiGenerateContent(prompt);
             catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                return;
            end
            disp(response);
    end
done = true;
end
