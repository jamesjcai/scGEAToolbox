function [done] = i_setllmmodel(~, ~)

preftagname = 'llmodelprovider';
if ~ispref('scgeatoolbox', preftagname)
    % answer = questdlg('LLM model has not been set up. Set it up?');
    % if ~strcmp(answer, 'Yes'), return; end
    % [done] = ix_setwdpath(pathdefult);
else
    s = getpref('scgeatoolbox', preftagname);
    answer1 = questdlg(sprintf('%s', s), ...
        'LLM Model', ...
        'Use this', 'Use another', 'Cancel', 'Use this');
    switch answer1
        case 'Use this'
            done = true;
            return;
        case 'Use another'
            % [done] = ix_setwdpath(s);
        otherwise
            return;
    end
end



listItems = {'Ollama', 'OpenAI', 'DeepSeek', 'xAI'};
[selectedIndex, ok] = listdlg('PromptString', 'Select a LLM provider:', ...
                              'SelectionMode', 'single', ...
                              'ListString', listItems, ...
                              'ListSize', [220 300], ...
                              'InitialValue', 1);
if ok
    selectedProvider = listItems{selectedIndex};
    % fprintf('You selected: %s\n', selectedProvider);
    switch selectedProvider
        case 'Ollama'
            a = webread("http://localhost:11434");
            if strcmp(a, 'Ollama is running')
                [a,str]=dos('Ollama list');
                if a == 0
                tokens = regexp(str, '([a-zA-Z0-9.-]+):latest', 'tokens');
                model_names = cellfun(@(x) x{1}, tokens, 'UniformOutput', false);
                

                [idx, ok2] = listdlg('PromptString', 'Select a model:', ...
                              'SelectionMode', 'single', ...
                              'ListString', model_names, ...
                              'ListSize', [220 300]);
                    if ok2
                        selectedModel = model_names{idx};                        
                        setpref('scgeatoolbox', preftagname, ...
                            selectedProvider+":"+selectedModel);
                        done = true;
                    end
                end
            else
                helpdlg('Ollama is not running.');
                return;
            end
        otherwise
            helpdlg('The function supporting API is under development.','');
            return;
    end
else    
    % fprintf('No selection made.\n');
    return;
end

if done
     waitfor(helpdlg("LLM provider and model are set successfully.", ''));
end
