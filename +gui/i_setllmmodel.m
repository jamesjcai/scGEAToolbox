function [done] = i_setllmmodel(src, ~)

if nargin<1, src = []; end

[parentfig, ~] = gui.gui_getfigsce(src);
done = false;
preftagname = 'llapikeyenvfile';
if ~ispref('scgeatoolbox', preftagname)
    if ~strcmp('Yes', gui.myQuestdlg(parentfig, 'Locate LLM API key env file?')), return; end
    [file, path] = uigetfile('llm_api_key.env', 'Select File');
    if isequal(file, 0), return; end
    apikeyfile = fullfile(path, file);
    if isfile(apikeyfile)
        setpref('scgeatoolbox', preftagname, apikeyfile);
        if ~strcmp('Yes', gui.myQuestdlg(parentfig, "LLM API key env file is located successfully. Continue?"))
            return;
        end
    else
        gui.myHelpdlg(parentfig, "Invalid file.")
        return;
    end
else
    apikeyfile = getpref('scgeatoolbox', preftagname);
    answer1 = gui.myQuestdlg(parentfig, sprintf('%s', apikeyfile), ...
        'Selected API Key File', ...
        {'Use this', 'Use another', '🌐Learn api_key file...'}, 'Use this');
    if isempty(answer1), return; end
    switch answer1
        case 'Cancel'
            return;
        case 'Use another'
            [file, path] = uigetfile('llm_api_key.env', 'Select API Key File');
            if isequal(file, 0), return; end
            apikeyfile = fullfile(path, file);
            setpref('scgeatoolbox', preftagname, apikeyfile);
        case '🌐Learn api_key file...'
            pause(1)
            web('https://github.com/jamesjcai/scGEAToolbox/blob/main/assets/Misc/.env.example')
            return;
    end
end

preftagname = 'llmodelprovider';
if ispref('scgeatoolbox', preftagname)
    s = getpref('scgeatoolbox', preftagname);
    answer1 = gui.myQuestdlg(parentfig, sprintf('%s', s), ...
        'Selected LLM Model', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    if isempty(answer1), return; end
    switch answer1
        case 'Use this'
            fw = gui.myWaitbar(parentfig);
            [done] = llm.i_checkllm(apikeyfile);
            gui.myWaitbar(parentfig, fw);
            if done
                gui.myHelpdlg(parentfig, "LLM provider and" + ...
                 " model are set successfully.");
            else
                gui.myWarndlg(parentfig, "LLM provider and" + ...
                 " model are not set successfully.");
            end
            return;
        case 'Use another'
            
        otherwise
            return;
    end
end

listItems = {'Ollama', 'Gemini', 'TAMUAIChat', 'OpenAI', 'Anthropic', ...
             'DeepSeek', 'xAI', 'Mistral', 'Cohere'};

if gui.i_isuifig(parentfig)
    [selectedIndex, ok] = gui.myListdlg(parentfig, listItems, ...
            'Select a LLM provider:', listItems(1));
else
    [selectedIndex, ok] = listdlg('PromptString', ...
                          'Select a LLM provider:', ...
                          'SelectionMode', 'single', ...
                          'ListString', listItems, ...
                          'ListSize', [220 300], ...
                          'InitialValue', 1);
end

if ~ok, return; end

selectedProvider = listItems{selectedIndex};
switch selectedProvider
    case 'Ollama'
        a = '';
        try 
            a = webread("http://localhost:11434");
        catch
        end
        if strcmp(a, 'Ollama is running')
            [a,str]=dos('Ollama list');
            if a == 0
            tokens = regexp(str, '([a-zA-Z0-9.-]+):latest', 'tokens');
            model_names = cellfun(@(x) x{1}, tokens, 'UniformOutput', false);
            
           if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select a model:');
            else
                [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                              'SelectionMode', 'single', ...
                              'ListString', model_names, ...
                              'ListSize', [220 300]);
            end

                if ok2
                    selectedModel = model_names{idx};                        
                    setpref('scgeatoolbox', preftagname, ...
                        selectedProvider+":"+selectedModel);
                    done = true;
                else
                    return;
                end
            end
        else
            gui.myHelpdlg(parentfig, 'Ollama is not running.');
            return;
        end
    case 'Gemini'
        if ~exist(apikeyfile,"file")
            gui.myErrordlg(parentfig,"llm_api_key.env is not a valid file.");
        end
        loadenv(apikeyfile,"FileType","env");
        if ~isempty(getenv("GEMINI_API_KEY"))        
            url = sprintf('https://generativelanguage.googleapis.com/v1beta/models?key=%s', ...
                  getenv("GEMINI_API_KEY"));
            a = webread(url);
            model_names = cellfun(@(x) x.name, a.models, 'UniformOutput', false);
            model_names = extractAfter(model_names, 7);
            [y, idx]=ismember('gemini-2.0-flash', model_names);
            if y
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', model_names(idx));
                else                    
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300], 'InitialValue', idx);
                end
            else
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:');
                else                
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300]);
                end
            end
            if ok2
                selectedModel = model_names{idx};                        
                setpref('scgeatoolbox', preftagname, ...
                    selectedProvider+":"+selectedModel);
                done = true;
            else
                return;
            end
        end
    case 'TAMUAIChat'
        if ~exist(apikeyfile,"file")
            gui.myErrordlg(parentfig,"llm_api_key.env is not a valid file.");
        end
        loadenv(apikeyfile,"FileType","env");
        if ~isempty(getenv("TAMUAI_API_KEY"))        
            OPEN_WEBUI_API_ENDPOINT = "https://chat-api.tamu.ai";
            models_url = sprintf('%s/api/models', OPEN_WEBUI_API_ENDPOINT);
            
            options = weboptions('HeaderFields', {'Authorization', ...
                sprintf('Bearer %s', getenv("TAMUAI_API_KEY"))}, ...
                'ContentType', 'json', 'Timeout', 50);

            fw = gui.myWaitbar(parentfig);
            
            try
                models_response = webread(models_url, options);
                model_names = string(cellfun(@(s) s.id, models_response.data, 'UniformOutput', false));
            catch ME
                gui.myWaitbar(parentfig, fw, true);
                gui.myErrordlg(parentfig, ME.message, 'Error fetching models');
                return;
            end

            gui.myWaitbar(parentfig, fw);

            [y, idx]=ismember('protected.gpt-4.1', model_names);
            if y
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', model_names(idx));
                else
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300], 'InitialValue', idx);
                end
            else
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:');
                else                
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300]);
                end
            end
            if ok2
                selectedModel = model_names{idx};                        
                setpref('scgeatoolbox', preftagname, ...
                    selectedProvider+":"+selectedModel);
                done = true;
            else
                return;
            end
        end
    case 'OpenAI'
        if ~exist(apikeyfile,"file")
            gui.myErrordlg(parentfig,"llm_api_key.env is not a valid file.");
        end
        loadenv(apikeyfile,"FileType","env");
        if ~isempty(getenv("OpenAI_API_KEY"))        
            OPEN_WEBUI_API_ENDPOINT = "https://api.openai.com/v1";
            models_url = sprintf('%s/models', OPEN_WEBUI_API_ENDPOINT);
            
            options = weboptions('HeaderFields', {'Authorization', ...
                sprintf('Bearer %s', getenv("OpenAI_API_KEY"))}, ...
                'ContentType', 'json',...
                'Timeout', 30);
            
            try
                models_response = webread(models_url, options);
                model_names = string(arrayfun(@(s) s.id, models_response.data, 'UniformOutput', false));
             catch ME
                fprintf('Error fetching models: %s\n', ME.message);
                return;
            end

            [y, idx]=ismember('gpt-4.1', model_names);
            if y
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', model_names(idx));
                else
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300], 'InitialValue', idx);
                end
            else
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:');
                else               
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300]);
                end
            end
            if ok2
                selectedModel = model_names{idx};                        
                setpref('scgeatoolbox', preftagname, ...
                    selectedProvider+":"+selectedModel);
                done = true;
            else
                return;
            end
        end

    case 'Anthropic'
        % Anthropic Claude models via the official Messages API
        % Requires ANTHROPIC_API_KEY to be set in the env file.
        if ~exist(apikeyfile, "file")
            gui.myErrordlg(parentfig, "llm_api_key.env is not a valid file.");
            return;
        end
        loadenv(apikeyfile, "FileType", "env");
        api_key = getenv("ANTHROPIC_API_KEY");
        if isempty(api_key)
            gui.myWarndlg(parentfig, ...
                "ANTHROPIC_API_KEY not found in the env file. " + ...
                "Please add it and try again.");
            return;
        end

        % Fetch available models from the Anthropic API
        models_url = 'https://api.anthropic.com/v1/models';
        options = weboptions( ...
            'HeaderFields', { ...
                'x-api-key',         api_key; ...
                'anthropic-version', '2023-06-01'}, ...
            'ContentType', 'json', ...
            'Timeout', 30);

        fw = gui.myWaitbar(parentfig);
        try
            models_response = webread(models_url, options);
            % Response shape: struct with field 'data', each element has 'id'
            model_names = string(cellfun(@(s) s.id, ...
                models_response.data, 'UniformOutput', false));
        catch ME
            gui.myWaitbar(parentfig, fw, true);
            gui.myErrordlg(parentfig, ME.message, 'Error fetching Anthropic models');
            return;
        end
        gui.myWaitbar(parentfig, fw);

        % Pre-select claude-sonnet-4-6 as the recommended default if present
        preferred_model = 'claude-sonnet-4-6';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y
            % Fall back to first model
            idx = 1;
            y   = ~isempty(model_names);
        end

        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a Claude model:', model_names(idx));
            else
                [idx, ok2] = listdlg('PromptString', 'Select a Claude model:', ...
                              'SelectionMode', 'single', ...
                              'ListString', model_names, ...
                              'ListSize', [300 300], ...
                              'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a Claude model:');
            else
                [idx, ok2] = listdlg('PromptString', 'Select a Claude model:', ...
                              'SelectionMode', 'single', ...
                              'ListString', model_names, ...
                              'ListSize', [300 300]);
            end
        end

        if ok2
            selectedModel = model_names{idx};
            setpref('scgeatoolbox', preftagname, ...
                selectedProvider + ":" + selectedModel);
            done = true;
        else
            return;
        end

    case 'DeepSeek'
        % DeepSeek — OpenAI-compatible API
        % Requires DEEPSEEK_API_KEY in the env file.
        if ~exist(apikeyfile, "file")
            gui.myErrordlg(parentfig, "llm_api_key.env is not a valid file.");
            return;
        end
        loadenv(apikeyfile, "FileType", "env");
        api_key = getenv("DEEPSEEK_API_KEY");
        if isempty(api_key)
            gui.myWarndlg(parentfig, ...
                "DEEPSEEK_API_KEY not found in the env file. " + ...
                "Please add it and try again.");
            return;
        end
        models_url = 'https://api.deepseek.com/v1/models';
        options = weboptions('HeaderFields', {'Authorization', ...
            sprintf('Bearer %s', api_key)}, ...
            'ContentType', 'json', 'Timeout', 30);
        fw = gui.myWaitbar(parentfig);
        try
            models_response = webread(models_url, options);
            model_names = string(cellfun(@(s) s.id, ...
                models_response.data, 'UniformOutput', false));
        catch ME
            gui.myWaitbar(parentfig, fw, true);
            gui.myErrordlg(parentfig, ME.message, 'Error fetching DeepSeek models');
            return;
        end
        gui.myWaitbar(parentfig, fw);
        preferred_model = 'deepseek-chat';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y, idx = 1; y = ~isempty(model_names); end
        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a DeepSeek model:', model_names(idx));
            else
                [idx, ok2] = listdlg('PromptString', 'Select a DeepSeek model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300], 'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select a DeepSeek model:');
            else
                [idx, ok2] = listdlg('PromptString', 'Select a DeepSeek model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300]);
            end
        end
        if ok2
            selectedModel = model_names{idx};
            setpref('scgeatoolbox', preftagname, selectedProvider + ":" + selectedModel);
            done = true;
        else
            return;
        end

    case 'xAI'
        % xAI Grok — OpenAI-compatible API
        % Requires XAI_API_KEY in the env file.
        if ~exist(apikeyfile, "file")
            gui.myErrordlg(parentfig, "llm_api_key.env is not a valid file.");
            return;
        end
        loadenv(apikeyfile, "FileType", "env");
        api_key = getenv("XAI_API_KEY");
        if isempty(api_key)
            gui.myWarndlg(parentfig, ...
                "XAI_API_KEY not found in the env file. " + ...
                "Please add it and try again.");
            return;
        end
        models_url = 'https://api.x.ai/v1/models';
        options = weboptions('HeaderFields', {'Authorization', ...
            sprintf('Bearer %s', api_key)}, ...
            'ContentType', 'json', 'Timeout', 30);
        fw = gui.myWaitbar(parentfig);
        try
            models_response = webread(models_url, options);
            model_names = string(cellfun(@(s) s.id, ...
                models_response.data, 'UniformOutput', false));
        catch ME
            gui.myWaitbar(parentfig, fw, true);
            gui.myErrordlg(parentfig, ME.message, 'Error fetching xAI models');
            return;
        end
        gui.myWaitbar(parentfig, fw);
        preferred_model = 'grok-3';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y, idx = 1; y = ~isempty(model_names); end
        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select an xAI model:', model_names(idx));
            else
                [idx, ok2] = listdlg('PromptString', 'Select an xAI model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300], 'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select an xAI model:');
            else
                [idx, ok2] = listdlg('PromptString', 'Select an xAI model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300]);
            end
        end
        if ok2
            selectedModel = model_names{idx};
            setpref('scgeatoolbox', preftagname, selectedProvider + ":" + selectedModel);
            done = true;
        else
            return;
        end

    case 'Mistral'
        % Mistral AI — OpenAI-compatible API
        % Requires MISTRAL_API_KEY in the env file.
        if ~exist(apikeyfile, "file")
            gui.myErrordlg(parentfig, "llm_api_key.env is not a valid file.");
            return;
        end
        loadenv(apikeyfile, "FileType", "env");
        api_key = getenv("MISTRAL_API_KEY");
        if isempty(api_key)
            gui.myWarndlg(parentfig, ...
                "MISTRAL_API_KEY not found in the env file. " + ...
                "Please add it and try again.");
            return;
        end
        models_url = 'https://api.mistral.ai/v1/models';
        options = weboptions('HeaderFields', {'Authorization', ...
            sprintf('Bearer %s', api_key)}, ...
            'ContentType', 'json', 'Timeout', 30);
        fw = gui.myWaitbar(parentfig);
        try
            models_response = webread(models_url, options);
            model_names = string(cellfun(@(s) s.id, ...
                models_response.data, 'UniformOutput', false));
            % Keep only chat-capable models (exclude embed/moderation models)
            is_chat = cellfun(@(s) isfield(s,'capabilities') && ...
                isfield(s.capabilities,'completion_chat') && ...
                s.capabilities.completion_chat, models_response.data);
            if any(is_chat)
                model_names = model_names(is_chat);
            end
        catch ME
            gui.myWaitbar(parentfig, fw, true);
            gui.myErrordlg(parentfig, ME.message, 'Error fetching Mistral models');
            return;
        end
        gui.myWaitbar(parentfig, fw);
        preferred_model = 'mistral-large-latest';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y, idx = 1; y = ~isempty(model_names); end
        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a Mistral model:', model_names(idx));
            else
                [idx, ok2] = listdlg('PromptString', 'Select a Mistral model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300], 'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select a Mistral model:');
            else
                [idx, ok2] = listdlg('PromptString', 'Select a Mistral model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300]);
            end
        end
        if ok2
            selectedModel = model_names{idx};
            setpref('scgeatoolbox', preftagname, selectedProvider + ":" + selectedModel);
            done = true;
        else
            return;
        end

    case 'Cohere'
        % Cohere — native API (slightly different from OpenAI-compatible)
        % Requires COHERE_API_KEY in the env file.
        if ~exist(apikeyfile, "file")
            gui.myErrordlg(parentfig, "llm_api_key.env is not a valid file.");
            return;
        end
        loadenv(apikeyfile, "FileType", "env");
        api_key = getenv("COHERE_API_KEY");
        if isempty(api_key)
            gui.myWarndlg(parentfig, ...
                "COHERE_API_KEY not found in the env file. " + ...
                "Please add it and try again.");
            return;
        end
        % Use /v2/models?endpoint=chat to get only chat-capable models
        models_url = 'https://api.cohere.com/v2/models?endpoint=chat&page_size=50';
        options = weboptions('HeaderFields', { ...
            'Authorization', sprintf('Bearer %s', api_key); ...
            'X-Client-Name', 'scGEAToolbox'}, ...
            'ContentType', 'json', 'Timeout', 30);
        fw = gui.myWaitbar(parentfig);
        try
            models_response = webread(models_url, options);
            % Cohere response: struct with field 'models', each has 'name'
            model_names = string(cellfun(@(s) s.name, ...
                models_response.models, 'UniformOutput', false));
        catch ME
            gui.myWaitbar(parentfig, fw, true);
            gui.myErrordlg(parentfig, ME.message, 'Error fetching Cohere models');
            return;
        end
        gui.myWaitbar(parentfig, fw);
        preferred_model = 'command-r-plus';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y, idx = 1; y = ~isempty(model_names); end
        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a Cohere model:', model_names(idx));
            else
                [idx, ok2] = listdlg('PromptString', 'Select a Cohere model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300], 'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select a Cohere model:');
            else
                [idx, ok2] = listdlg('PromptString', 'Select a Cohere model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300]);
            end
        end
        if ok2
            selectedModel = model_names{idx};
            setpref('scgeatoolbox', preftagname, selectedProvider + ":" + selectedModel);
            done = true;
        else
            return;
        end

    otherwise
        gui.myWarndlg(parentfig, ...
            sprintf(['The function supporting %s API is ' ...
            'under development.'], ...
            selectedProvider));
        return;
end

fw = gui.myWaitbar(parentfig);
% [done2] = llm.i_checkllm(apikeyfile);
done2=true;
gui.myWaitbar(parentfig, fw);

if done && done2
     gui.myHelpdlg(parentfig, "LLM provider and" + ...
         " model are set successfully.");
end
