function [done] = i_checkllm(apikeyfile, provider, parentfig)
done = false;

if nargin<3, parentfig = []; end
if nargin<2, provider = 'Gemini'; end
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

    preftagname = 'llmodelprovider';
    s = getpref('scgeatoolbox', preftagname);
    providermodel = strsplit(s,':');

    provider = providermodel{1};
    switch provider
        case 'Ollama'
            try
                chat = ollamaChat(providermodel{2}, TimeOut = 1200);
                prompt = "Why is the sky blue?";
                feedbk = generate(chat, prompt);
                done = true;
            catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
            end
        case 'TAMUAIChat'
            % Test a chat completion to llama3.2 model
            loadenv(apikeyfile,"FileType","env");
            OPEN_WEBUI_API_ENDPOINT = "https://chat-api.preview.tamu.ai";
            OPEN_WEBUI_API_KEY = getenv("TAMUAI_API_KEY");
            chat_url = sprintf('%s/api/chat/completions', OPEN_WEBUI_API_ENDPOINT);
            
            % Create request body structure
            messages_cell = {struct('role', 'user', 'content', 'Why is the sky blue?')};
            body_struct = struct('model', providermodel{2}, ...
                                'stream', false, ...
                                'messages', {messages_cell});            
            % Debug: Display the request body structure
            fprintf('Request body structure:\n');
            disp(body_struct);
            
            % Set up options for webwrite (POST request)
            post_options = weboptions('MediaType', 'application/json', ...
                                     'RequestMethod', 'POST', ...
                                     'HeaderFields', {'Authorization', sprintf('Bearer %s', OPEN_WEBUI_API_KEY)});
            
            % Make the chat completion request
            try
                chat_response = webwrite(chat_url, body_struct, post_options);
                % Display chat response as JSON
                chat_json = jsonencode(chat_response);
                fprintf('Chat completion response:\n%s\n', chat_json);
                done = true;
            catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                % Additional debug info
                if contains(ME.message, '400')
                    fprintf('This is likely a request format issue. Check the API documentation for the exact expected format.\n');
                end
            end   
        case 'Gemini'
            loadenv(apikeyfile,"FileType","env");
            apiKey = getenv("GEMINI_API_KEY");
            model = providermodel{2};

            prompt = "Why is the sky blue?";
            try
                response = llm.callGemini2(apiKey, prompt, model);
                % response = llm.geminiGenerateContent(prompt);
             catch ME
                fprintf('Error in chat completion: %s\n', ME.message);
                if contains(ME.message, '400')
                    fprintf('This is likely a request format issue. Check the API documentation for the exact expected format.\n');
                end
            end   

    end
