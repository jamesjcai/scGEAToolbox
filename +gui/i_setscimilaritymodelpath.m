function [selectedDir] = i_setscimilaritymodelpath

selectedDir = '';
preftagname = 'scimilmodelpath';
if ispref('scgeatoolbox', preftagname)
    selectedDir = getpref('scgeatoolbox', preftagname);
end

if isempty(selectedDir) || ~isfolder(selectedDir)
    answer = questdlg('Scimilarity model path has not been set up. Locate it?');
    if strcmp('Yes', answer)
        [done] = ix_setpath;
        if ~done, return; end
        waitfor(helpdlg('Scimilarity model path is set successfully.', ''));
    else
        return;
    end
else
    answer = questdlg(sprintf('%s', selectedDir), ...
        'Model Path', ...
        'Use this', 'Use another', 'Cancel', 'Use this');
    switch answer
        case 'Use this'
            done = true;
        case 'Use another'
            if ~ix_setpath
                return;
            end
            done = true;
            waitfor(helpdlg('Scimilarity model path is set successfully.', ''));
        case {'Cancel', ''}
            selectedDir = '';
            done = false;
        otherwise
            selectedDir = '';
            done = false;
    end
end

if ~done
    warndlg('SCimilarity model path is not set.', '');
end


    function [done] = ix_setpath
        done = false;
        promptTitle = 'Select a folder that contains the model';
        selectedDir = uigetdir(pwd, promptTitle);
        if selectedDir == 0
            fprintf('Folder selection canceled.\n');
            selectedDir = '';            
        else
            fprintf('Selected folder: %s\n', selectedDir);
            done = true;
            setpref('scgeatoolbox', preftagname, selectedDir);
        end
    end

end