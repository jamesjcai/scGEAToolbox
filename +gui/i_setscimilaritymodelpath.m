function [selectedDir] = i_setscimilaritymodelpath(src, ~, parentfig)

if nargin<3, parentfig = []; end
if ~isempty(src)
    [parentfig] = gui.gui_getfigsce(src);
end
selectedDir = '';
preftagname = 'scimilmodelpath';
if ispref('scgeatoolbox', preftagname)
    selectedDir = getpref('scgeatoolbox', preftagname);
end

if isempty(selectedDir) || ~isfolder(selectedDir)
    answer = gui.myQuestdlg(parentfig, 'Scimilarity model path has not been set up. Locate it?');
    if strcmp('Yes', answer)
        [done] = ix_setpath;
        if ~done, return; end
        gui.myHelpdlg(parentfig, 'Scimilarity model path is set successfully.');
    else
        return;
    end
else
    answer = gui.myQuestdlg(parentfig, sprintf('%s', selectedDir), ...
        'Model Path', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    switch answer
        case 'Use this'
            %done = true;
        case 'Use another'
            if ~ix_setpath
                return;
            end
            %done = true;
            gui.myHelpdlg(parentfig, ...
                'Scimilarity model path is set successfully.');
            return;
        case {'Cancel', ''}
            selectedDir = '';
            %done = false;
        otherwise
            selectedDir = '';
            %done = false;
    end
end

%if ~done && (isempty(selectedDir) || ~isfolder(selectedDir))
%    gui.myWarndlg(parentfig, 'SCimilarity model path is not set.');
%end


    function [y] = ix_setpath
        y = false;
        promptTitle = 'Select a folder that contains the model';
        selectedDir = uigetdir(pwd, promptTitle);
        if isvalid(parentfig)
            figure(parentfig);
        end
        if selectedDir == 0
            fprintf('Folder selection canceled.\n');
            selectedDir = '';            
        else
            fprintf('Selected folder: %s\n', selectedDir);
            y = true;
            setpref('scgeatoolbox', preftagname, selectedDir);
        end

        label_ints_file = fullfile(selectedDir, 'label_ints.csv');
        if ~exist(label_ints_file, "file")
            y = false;
        end
    end

end