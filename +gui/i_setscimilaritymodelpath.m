function [selectedDir] = i_setscimilaritymodelpath(src, ~, parentfig)

if nargin<3, parentfig = []; end
if ~isempty(src)
    [parentfig] = gui.gui_getfigsce(src);
end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
selectedDir = '';
preftagname = 'scimilmodelpath';
if ispref('scgeatoolbox', preftagname)
    selectedDir = getpref('scgeatoolbox', preftagname);
end

if isempty(selectedDir) || ~isfolder(selectedDir)

    % Locate is offered straight away. It used to sit behind a Download
    % question whose No ended the function, so a user who already had
    % the model had to open Zenodo to reach the folder picker.
    answer = gui.myQuestdlg(parentfig, ['The SCimilarity model folder ' ...
        'is not set up. Locate a downloaded model, or download it first ' ...
        '(a large tarball; downloading and unpacking take several minutes)?'], ...
        'Model Path', {'Locate Folder', 'Download', 'Cancel'}, 'Locate Folder');
    switch answer
        case 'Locate Folder'
        case 'Download'
            web('https://zenodo.org/records/10685499');
            if ~strcmp('Yes', gui.myQuestdlg(parentfig, ...
                    'Once the model is downloaded and unpacked, locate its folder?'))
                return;
            end
        otherwise
            return;
    end
    if ~ix_setpath, return; end
    gui.myHelpdlg(parentfig, 'Scimilarity model path is set successfully.');
else
    answer = gui.myQuestdlg(parentfig, sprintf('%s', selectedDir), ...
        'Model Path', ...
        {'Use this', 'Use another', 'Cancel'}, 'Use this');
    switch answer
        case 'Use this'
        case 'Use another'
            if ~ix_setpath
                return;
            end
            gui.myHelpdlg(parentfig, ...
                'Scimilarity model path is set successfully.');
            return;
        case {'Cancel', ''}
            selectedDir = '';
        otherwise
            selectedDir = '';
    end
end

% if ~done && (isempty(selectedDir) || ~isfolder(selectedDir))
%    gui.myWarndlg(parentfig, 'SCimilarity model path is not set.');
% end


function [y] = ix_setpath
        y = false;
        promptTitle = 'Select a folder that contains the model';
        selectedDir = uigetdir(pwd, promptTitle);
        if pkg.i_isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure')
            figure(parentfig);
        end

        if isequal(selectedDir, 0)
            fprintf('Folder selection canceled.\n');
            selectedDir = '';
            return;
        end
        % Checked before it is saved. The folder used to be stored first,
        % so a wrong one was kept and offered as "Use this" next time.
        if ~isfile(fullfile(selectedDir, 'label_ints.csv'))
            gui.myWarndlg(parentfig, sprintf(['%s does not look like an ' ...
                'SCimilarity model folder (no label_ints.csv). The path ' ...
                'was not saved.'], selectedDir));
            selectedDir = '';
            return;
        end
        fprintf('Selected folder: %s\n', selectedDir);
        setpref('scgeatoolbox', preftagname, selectedDir);
        y = true;
    end

end
