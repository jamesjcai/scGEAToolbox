function [done] = i_setllmmodel(src, ~)

if nargin<1, src = []; end

[parentfig, ~] = gui.gui_getfigsce(src);
done = false;
% Keys come from the environment, the MATLAB vault or an llm_api_key.env
% file (llm.i_getapikey). The env file is optional: choosing the vault
% stores an empty preference, and each provider below asks for its key,
% through MATLAB's masked SETSECRET dialog, the first time it is missing.
preftagname = 'llapikeyenvfile';
apikeyfile = "";
if ~ispref('scgeatoolbox', preftagname)
    answer0 = gui.myQuestdlg(parentfig, ...
        ['Keep LLM API keys encrypted in the MATLAB vault, ' ...
        'or read them from an llm_api_key.env file?'], ...
        'LLM API Keys', {'MATLAB vault', 'Env file...', 'Cancel'}, 'MATLAB vault');
    switch answer0
        case 'MATLAB vault'
            setpref('scgeatoolbox', preftagname, '');
        case 'Env file...'
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
        otherwise
            return;
    end
elseif isempty(getpref('scgeatoolbox', preftagname))
    % The vault was chosen before: no file to confirm.
else
    apikeyfile = getpref('scgeatoolbox', preftagname);
    answer1 = gui.myQuestdlg(parentfig, sprintf('%s', apikeyfile), ...
        'Selected API Key File', ...
        {'Use this', 'Use another', '🌐Learn api_key file...'}, 'Use this');
    if isempty(answer1), return; end
    switch answer1
        case 'Use another'
            [file, path] = uigetfile('llm_api_key.env', 'Select API Key File');
            if isequal(file, 0), return; end
            apikeyfile = fullfile(path, file);
            setpref('scgeatoolbox', preftagname, apikeyfile);
        case '🌐Learn api_key file...'
            drawnow;
            web('https://github.com/jamesjcai/scGEAToolbox/blob/main/assets/Misc/.env.example')
            return;
    end
end
fileValues = llm.i_readkeyfile(apikeyfile);

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
            [done, tested] = llm.i_checkllm(apikeyfile);
            gui.myWaitbar(parentfig, fw);
            if done && ~tested
                % Said so rather than reported as a success: the check has
                % no test call for every provider.
                gui.myHelpdlg(parentfig, "LLM provider and model are " + ...
                    "kept, but this provider cannot be tested from here.");
            elseif done
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

% Only providers some feature can call. Anthropic and Cohere used to be
% offered too, but neither llm.i_askllm (the LLM reports) nor GEOcellar has
% a client for them, so choosing one broke both. The two features support
% different sets, and each entry says which it works with.
listItems = {'Ollama', 'Gemini', 'TAMUAIChat', 'OpenAI', ...
             'NVIDIA', 'DeepSeek', 'xAI', 'Mistral'};
listShown = {'Ollama', 'Gemini (LLM reports only, not GEOcellar)', ...
             'TAMUAIChat', 'OpenAI', 'NVIDIA', ...
             'DeepSeek (GEOcellar only, not LLM reports)', ...
             'xAI (GEOcellar only, not LLM reports)', ...
             'Mistral (GEOcellar only, not LLM reports)'};

if gui.i_isuifig(parentfig)
    [selectedIndex, ok] = gui.myListdlg(parentfig, listShown, ...
            'Select a LLM provider:', listShown(1), false);
else
    [selectedIndex, ok] = listdlg('PromptString', ...
                          'Select a LLM provider:', ...
                          'SelectionMode', 'single', ...
                          'ListString', listShown, ...
                          'ListSize', [220 300], ...
                          'InitialValue', 1);
end

if ~ok, return; end

selectedProvider = listItems{selectedIndex};
switch selectedProvider
    case 'Ollama'
        a = '';
        try
            % Use 127.0.0.1 rather than localhost: Ollama binds 127.0.0.1, but
            % localhost resolves to ::1 first on some hosts, which would make a
            % working install look absent.
            a = webread("http://127.0.0.1:11434");
        catch
            % ollama may not be running; a stays empty and the check below fails over
        end
        if strcmp(a, 'Ollama is running')
            a = webread('http://127.0.0.1:11434/api/tags');


           if isfield(a, 'models') && size(a.models, 1) > 0
               % jsondecode returns a struct array when every model object has
               % the same fields and a cell array when they differ, so both
               % shapes have to be handled.
               if iscell(a.models)
                   model_names = string(cellfun(@(m) m.name, a.models, ...
                       'UniformOutput', false));
               else
                   model_names = string({a.models.name});
               end

           if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                    'Select a model:', [], false);
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
            else
                gui.myWarndlg(parentfig, sprintf([ ...
                    'Ollama is running, but no models are installed.\n\n' ...
                    'Pull one from a terminal first, for example:\n' ...
                    '    ollama pull qwen3:4b']));
                return;
            end
        else
            gui.myHelpdlg(parentfig, 'Ollama is not running.');
            return;
        end
    case 'Gemini'
        apiKey = i_providerkey(parentfig, "GEMINI_API_KEY", fileValues);
        if ~isempty(apiKey)
            % The key in a header, not in the URL, where it is written to
            % proxy and server access logs.
            % Gemini answers a bad key with 400 API_KEY_INVALID, not 401.
            url = 'https://generativelanguage.googleapis.com/v1beta/models';
            [a, ok] = i_fetchmodels(parentfig, url, ...
                @(k) weboptions('HeaderFields', {'x-goog-api-key', k}, 'Timeout', 30), ...
                "GEMINI_API_KEY", apiKey, fileValues, "Gemini", [400 401 403]);
            if ~ok, return; end
            model_names = cellstr(i_field2str(a.models, 'name'));
            model_names = extractAfter(model_names, 7);
            [y, idx]=ismember('gemini-2.0-flash', model_names);
            if y
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', model_names(idx), false);
                else
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300], 'InitialValue', idx);
                end
            else
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', [], false);
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
    case 'NVIDIA'
        apiKey = i_providerkey(parentfig, "NVIDIA_API_KEY", fileValues);
        if ~isempty(apiKey)
            OPEN_WEBUI_API_ENDPOINT = i_setting("NVIDIA_API_BASE", ...
                "https://integrate.api.nvidia.com/v1", fileValues);
            models_url = sprintf('%s/models', OPEN_WEBUI_API_ENDPOINT);

            [models_response, ok] = i_fetchmodels(parentfig, models_url, ...
                @(k) i_beareroptions(k, 50), "NVIDIA_API_KEY", apiKey, fileValues, "NVIDIA");
            if ~ok, return; end
            model_names = i_field2str(models_response.data, 'id');

            [y, idx]=ismember('minimaxai/minimax-m2.5', model_names);
            if y
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a model:', model_names(idx), false);
                else
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                        'SelectionMode', 'single', ...
                        'ListString', model_names, ...
                        'ListSize', [220 300], 'InitialValue', idx);
                end
            else
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a model:', [], false);
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
        apiKey = i_providerkey(parentfig, "TAMUAI_API_KEY", fileValues);
        if ~isempty(apiKey)
            OPEN_WEBUI_API_ENDPOINT = "https://chat-api.tamu.ai";
            models_url = sprintf('%s/api/models', OPEN_WEBUI_API_ENDPOINT);

            [models_response, ok] = i_fetchmodels(parentfig, models_url, ...
                @(k) i_beareroptions(k, 50), "TAMUAI_API_KEY", apiKey, fileValues, "TAMU AI");
            if ~ok, return; end
            model_names = i_field2str(models_response.data, 'id');

            [y, idx]=ismember('protected.gpt-4.1', model_names);  % xxx
            if y
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', model_names(idx), false);
                else
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300], 'InitialValue', idx);
                end
            else
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', [], false);
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
        % OPENAI_API_KEY, the name the pipelines and .env.example use, or
        % the OpenAI_API_KEY this once read alone: environment names are
        % case-sensitive on macOS and Linux, so there the conventional name
        % was never found and choosing OpenAI silently did nothing.
        openaiKey = i_providerkey(parentfig, ["OPENAI_API_KEY", "OpenAI_API_KEY"], fileValues);
        if ~isempty(openaiKey)
            OPEN_WEBUI_API_ENDPOINT = "https://api.openai.com/v1";
            models_url = sprintf('%s/models', OPEN_WEBUI_API_ENDPOINT);

            [models_response, ok] = i_fetchmodels(parentfig, models_url, ...
                @(k) i_beareroptions(k, 30), ["OPENAI_API_KEY", "OpenAI_API_KEY"], ...
                openaiKey, fileValues, "OpenAI");
            if ~ok, return; end
            model_names = i_field2str(models_response.data, 'id');

            [y, idx]=ismember('gpt-4.1', model_names);
            if y
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', model_names(idx), false);
                else
                    [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                                  'SelectionMode', 'single', ...
                                  'ListString', model_names, ...
                                  'ListSize', [220 300], 'InitialValue', idx);
                end
            else
                if gui.i_isuifig(parentfig)
                    [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                            'Select a model:', [], false);
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

    case 'DeepSeek'
        % DeepSeek — OpenAI-compatible API
        % Requires DEEPSEEK_API_KEY in the env file or the MATLAB vault.
        api_key = i_providerkey(parentfig, "DEEPSEEK_API_KEY", fileValues);
        if isempty(api_key), return; end
        models_url = 'https://api.deepseek.com/v1/models';
        [models_response, ok] = i_fetchmodels(parentfig, models_url, ...
            @(k) i_beareroptions(k, 30), "DEEPSEEK_API_KEY", api_key, fileValues, "DeepSeek");
        if ~ok, return; end
        model_names = i_field2str(models_response.data, 'id');
        preferred_model = 'deepseek-chat';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y, idx = 1; y = ~isempty(model_names); end
        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a DeepSeek model:', model_names(idx), false);
            else
                [idx, ok2] = listdlg('PromptString', 'Select a DeepSeek model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300], 'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select a DeepSeek model:', [], false);
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
        % Requires XAI_API_KEY in the env file or the MATLAB vault.
        api_key = i_providerkey(parentfig, "XAI_API_KEY", fileValues);
        if isempty(api_key), return; end
        models_url = 'https://api.x.ai/v1/models';
        [models_response, ok] = i_fetchmodels(parentfig, models_url, ...
            @(k) i_beareroptions(k, 30), "XAI_API_KEY", api_key, fileValues, "xAI");
        if ~ok, return; end
        model_names = i_field2str(models_response.data, 'id');
        preferred_model = 'grok-3';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y, idx = 1; y = ~isempty(model_names); end
        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select an xAI model:', model_names(idx), false);
            else
                [idx, ok2] = listdlg('PromptString', 'Select an xAI model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300], 'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select an xAI model:', [], false);
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
        % Requires MISTRAL_API_KEY in the env file or the MATLAB vault.
        api_key = i_providerkey(parentfig, "MISTRAL_API_KEY", fileValues);
        if isempty(api_key), return; end
        models_url = 'https://api.mistral.ai/v1/models';
        [models_response, ok] = i_fetchmodels(parentfig, models_url, ...
            @(k) i_beareroptions(k, 30), "MISTRAL_API_KEY", api_key, fileValues, "Mistral");
        if ~ok, return; end
        model_names = i_field2str(models_response.data, 'id');
        % Keep only chat-capable models (exclude embed/moderation models).
        % NUM2CELL lets one CELLFUN cover both JSONDECODE shapes.
        data = models_response.data;
        if ~iscell(data), data = num2cell(data); end
        is_chat = cellfun(@(s) isfield(s,'capabilities') && ...
            isfield(s.capabilities,'completion_chat') && ...
            isequal(s.capabilities.completion_chat, true), data);
        if any(is_chat)
            model_names = model_names(is_chat);
        end
        preferred_model = 'mistral-large-latest';
        [y, idx] = ismember(preferred_model, model_names);
        if ~y, idx = 1; y = ~isempty(model_names); end
        if y
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, ...
                        'Select a Mistral model:', model_names(idx), false);
            else
                [idx, ok2] = listdlg('PromptString', 'Select a Mistral model:', ...
                              'SelectionMode', 'single', 'ListString', model_names, ...
                              'ListSize', [300 300], 'InitialValue', idx);
            end
        else
            if gui.i_isuifig(parentfig)
                [idx, ok2] = gui.myListdlg(parentfig, model_names, 'Select a Mistral model:', [], false);
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

    otherwise
        gui.myWarndlg(parentfig, ...
            sprintf(['The function supporting %s API is ' ...
            'under development.'], ...
            selectedProvider));
        return;
end

if done
     gui.myHelpdlg(parentfig, "LLM provider and" + ...
         " model are set successfully.");
end
end

function key = i_providerkey(parentfig, names, fileValues)
% The provider's API key (llm.i_getapikey). When there is none, offer to
% enter it in MATLAB's masked dialog and keep it in the vault; '' when the
% user declines. The branches used to stop at "not a valid file" whenever
% no env file was set, even with the key in the environment.
key = llm.i_getapikey(names, fileValues);
if ~isempty(key), return; end
name = string(names(1));
answer = gui.myQuestdlg(parentfig, sprintf(['%s was not found in the ' ...
    'environment, the MATLAB vault or the env file. Enter it now? It is ' ...
    'stored encrypted in the MATLAB vault.'], name));
if strcmp(answer, 'Yes')
    key = llm.i_storeapikey(name);
end
end

function [response, ok] = i_fetchmodels(parentfig, url, makeOptions, names, key, fileValues, label, authCodes)
% WEBREAD(URL, MAKEOPTIONS(KEY)) behind a waitbar. When the provider
% rejects the key (an HTTP status in AUTHCODES, 401 and 403 by default),
% offer to enter a new one and try again. A wrong key used to end at an
% error dialog, and since it stayed in the vault every later attempt
% failed the same way. Any other failure is reported and ends the attempt.
if nargin < 8, authCodes = [401 403]; end
response = [];
ok = false;
while true
    fw = gui.myWaitbar(parentfig);
    try
        response = webread(url, makeOptions(key));
        gui.myWaitbar(parentfig, fw);
        ok = true;
        return;
    catch ME
        gui.myWaitbar(parentfig, fw, true);
        status = str2double(regexp(ME.identifier, 'HTTP(\d+)StatusCodeError', ...
            'tokens', 'once'));
        if ~any(status == authCodes)
            gui.myErrordlg(parentfig, ME.message, "Error fetching " + label + " models");
            return;
        end
    end
    key = i_replacekey(parentfig, names, fileValues, label, status);
    if isempty(key), return; end
end
end

function key = i_replacekey(parentfig, names, fileValues, label, status)
% Ask for a replacement for a rejected key and store it in the vault; ''
% when declined. A key in the environment cannot be replaced from here:
% it comes before the vault, so the new key would never be used.
key = '';
name = string(names(1));
[~, source] = llm.i_getapikey(names, fileValues);
if source == "environment"
    gui.myErrordlg(parentfig, sprintf(['%s rejected the key (HTTP %d). ' ...
        'It is set as the environment variable %s, which comes before ' ...
        'the MATLAB vault. Correct or clear it there, for example with ' ...
        'setenv("%s", ""), and try again.'], label, status, name, name), ...
        'API Key Rejected');
    return;
end
if source == "file"
    where = 'It is stored in the MATLAB vault and used instead of the one in the env file.';
else
    where = 'It replaces the key in the MATLAB vault.';
end
answer = gui.myQuestdlg(parentfig, sprintf(['%s rejected the %s key ' ...
    '(HTTP %d). Enter a new key? %s'], label, name, status, where), ...
    'API Key Rejected');
if strcmp(answer, 'Yes')
    key = llm.i_storeapikey(name);
end
end

function options = i_beareroptions(key, timeout)
% WEBOPTIONS for the providers that take the key as a Bearer token.
options = weboptions('HeaderFields', {'Authorization', ...
    sprintf('Bearer %s', key)}, 'ContentType', 'json', 'Timeout', timeout);
end

function v = i_setting(name, default, fileValues)
% Setting NAME from the environment, else the env file, else DEFAULT.
v = getenv(name);
if isempty(v) && isConfigured(fileValues) && isKey(fileValues, name)
    v = char(fileValues(name));
end
if isempty(v)
    v = default;
end
end

function v = i_field2str(d, name)
% D.(NAME) of every element as a string array, whichever shape JSONDECODE
% gave D: a struct array when every object in the reply has the same fields,
% a cell array of structs when they differ. The model-list code assumed one
% shape per provider -- CELLFUN for most, ARRAYFUN for OpenAI -- and threw
% "Error fetching models" whenever the reply came in the other.
if iscell(d)
    v = string(cellfun(@(x) x.(name), d, 'UniformOutput', false));
else
    v = string({d.(name)});
end
v = v(:);
end
