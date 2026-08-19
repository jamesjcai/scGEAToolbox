function [versionStr] = i_get_versionnum

versionStr = '';
mfolder = fileparts(mfilename('fullpath'));
% mfolder is the absolute path to +pkg; its parent is the toolbox root.
rootFolder = fileparts(mfolder);
xfilelocal = fullfile(rootFolder, "doc", "helptoc.xml");

if isfile(xfilelocal)
    txt = fileread(xfilelocal);
    % Look for version in lines containing "v25.7.7" or similar
    tokens = regexp(txt, 'v?(\d+\.\d+\.\d+)', 'tokens');

    if ~isempty(tokens)
        versionStr = tokens{1}{1};
        % fprintf('Extracted version: %s\n', versionStr);
    else
        warning('Version string not found in the file.');
    end
end
end
