function [done] = i_setllmmodel(src, ~)

[parentfig, ~] = gui.gui_getfigsce(src);
done = false;
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
    answer1 = gui.myQuestdlg(parentfig, sprintf('%s', apikeyfile), ...
        'Selected API Key File', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    if isempty(answer1), return; end
    switch answer1
        case 'Cancel'
            return;
        case 'Use another'
            [file, path] = uigetfile('llm_api_key.env', 'Select API Key File');
            if isequal(file, 0), return; end
            apikeyfile = fullfile(path, file);
            setpref('scgeatoolbox', preftagname, apikeyfile);
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

% listItems = {'Ollama', 'Gemini', 'TAMUAIChat', 'OpenAI', 'DeepSeek', 'xAI'};
listItems = {'Ollama', 'Gemini', 'TAMUAIChat'};

if gui.i_isuifig(parentfig)
    [selectedIndex, ok] = gui.myListdlg(parentfig, listItems, ...
            'Select a LLM provider:', listItems(1));
else
    [selectedIndex, ok] = listdlg('PromptString', 'Select a LLM provider:', ...
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
            % assignin("base", "a", a);
            % gui.myHelpdlg(parentfig, "GEMINI_API_KEY is load successfully.");
            % model_names={'gemini-2.0-flash'};
            model_names = cellfun(@(x) x.name, a.models, 'UniformOutput', false);
            model_names = extractAfter(model_names, 7);
            [y, idx]=ismember('gemini-2.0-flash', model_names);
            if y
                [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                              'SelectionMode', 'single', ...
                              'ListString', model_names, ...
                              'ListSize', [220 300], 'InitialValue', idx);
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
            end
        end
    case 'TAMUAIChat'
        if ~exist(apikeyfile,"file")
            gui.myErrordlg(parentfig,"llm_api_key.env is not a valid file.");
        end
        loadenv(apikeyfile,"FileType","env");
        if ~isempty(getenv("TAMUAI_API_KEY"))        
            % Show models available in json format
            OPEN_WEBUI_API_ENDPOINT = "https://chat-api.preview.tamu.ai";
            models_url = sprintf('%s/api/models', OPEN_WEBUI_API_ENDPOINT);
            
            % Set up options for webread
            options = weboptions('HeaderFields', {'Authorization', sprintf('Bearer %s', getenv("TAMUAI_API_KEY"))});
            
            % Make the GET request
            try
                models_response = webread(models_url, options);
                % Display models response as JSON
                % assignin("base",'models_response',  models_response)
                model_names = string({models_response.data.id});

            catch ME
                fprintf('Error fetching models: %s\n', ME.message);
                return;
            end

            [y, idx]=ismember('protected.gpt-4.1', model_names);
            if y
                [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                              'SelectionMode', 'single', ...
                              'ListString', model_names, ...
                              'ListSize', [220 300], 'InitialValue', idx);
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
            end
        end

    otherwise
        gui.myWarndlg(parentfig, ...
            sprintf(['The function supporting %s API is ' ...
            'under development.'], ...
            selectedProvider));
        return;
end

fw = gui.myWaitbar(parentfig);
[done2] = llm.i_checkllm(apikeyfile);
gui.myWaitbar(parentfig, fw);

if done && done2
     gui.myHelpdlg(parentfig, "LLM provider and" + ...
         " model are set successfully.");
end
