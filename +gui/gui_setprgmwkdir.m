function [wrkdir] = gui_setprgmwkdir(extprogname, preftagname, parentfig)
if nargin<3, parentfig = []; end
wrkdir = '';
% extprogname = 'R_monocle3';
% preftagname = 'externalwrkpath';
if ~gui.i_setwrkdir(preftagname, parentfig), return; end
s = getpref('scgeatoolbox', preftagname, []);
if isempty(s)
    % GUI.I_SETWRKDIR said it was set up and it is not. Every caller here
    % tests `if isempty(wkdir), return; end`, so hand back empty rather
    % than throwing: an ERROR out of a menu callback reaches the user as an
    % unhandled red stack trace with nothing to act on.
    gui.myErrordlg(parentfig, ['The working folder preference is empty. ' ...
        'Set it with Setup > Set Working Folder, then try again.'], ...
        'gui:gui_setprgmwkdir:noWorkingPath');
    return;
end
s1 = sprintf('%s_workingfolder', extprogname);
wrkdir = fullfile(s, s1);

if ~exist(wrkdir,"dir")
    mkdir(wrkdir);
else
    answer = gui.myQuestdlg(parentfig, ...
        sprintf('%s existing. Overwrite?', wrkdir));
    if ~strcmp(answer,'Yes')
        wrkdir = '';
        return;
    else
        deleteAllFiles(wrkdir);
        % if ~strcmp('Yes', gui.myQuestdlg(parentfig, ...
        %        'Existing files in the working folder will be overwritten or deleted. Continue?'))
        %    wrkdir = '';
        %    return;
        % end
    end
end
fprintf('CURRENTWDIR = "%s"\n', wrkdir);


end


function deleteAllFiles(wkdir)
% Check if directory exists
if ~isfolder(wkdir)
    warning('Directory does not exist: %s', wkdir);
    return;
end

% Get all files in the directory (excluding subdirectories)
files = dir(fullfile(wkdir, '*'));
files = files(~[files.isdir]); % Remove directories from the list

% Check if there are any files
if ~isempty(files)
    fprintf('Found %s in %s\n', pkg.i_plural(length(files), 'file'), wkdir);

    % Delete all files
    for i = 1:length(files)
        filePath = fullfile(wkdir, files(i).name);
        try
            delete(filePath);
            fprintf('Deleted: %s\n', files(i).name);
        catch ME
            warning('Could not delete %s: %s', files(i).name, ME.message);
        end
    end
else
    fprintf('No files found in %s\n', wkdir);
end
end
