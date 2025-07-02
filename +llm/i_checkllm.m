function [done2] = i_checkllm(apikeyfile)
done2 = false;

    preftagname = 'llmodelprovider';
    s = getpref('scgeatoolbox', preftagname);
    providermodel = strsplit(s,':');
    switch providermodel{1}
        case 'Ollama'
            try
                chat = ollamaChat(providermodel{2}, TimeOut = 1200);
                prompt = "Why is the sky blue?";
                feedbk = generate(chat, prompt);
                done2 = true;
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
                done2 = true;
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
            % API endpoint
            url = ['https://generativelanguage.googleapis.com/v1beta/models/' model ':generateContent?key=' apiKey];
            prompt = "Why is the sky blue?";
            requestBody = struct('contents', struct('parts', struct('text', prompt)));
            jsonBody = jsonencode(requestBody);
            % Set up the HTTP request
            options = weboptions('ContentType', 'json', 'RequestMethod', 'post');
            % Make the API call
            try
                response = webwrite(url, jsonBody, options);
                % Extract the actual response text
                if isfield(response, 'candidates') && ~isempty(response.candidates) && ...
                   isfield(response.candidates{1}, 'content') && ...
                   isfield(response.candidates{1}.content, 'parts') && ...
                   ~isempty(response.candidates{1}.content.parts)
                    response = response.candidates{1}.content.parts{1}.text;
                    fprintf('Chat completion response:\n%s\n', response);
                    done2 = true;
                end
            catch e
                fprintf('Error calling Gemini API: %s\n', e.message);
            end

    end
