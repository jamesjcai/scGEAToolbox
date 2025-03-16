function [matchingFiles] = e_openallfiles(k1, k2)

% if nargin<1, k1='parentfig'; k2='myErrordlg'; end
if nargin<1, k1='parentfig'; k2='fw = gui.myWaitbar'; end


mfolder = fileparts(mfilename('fullpath'));
% Define the search directory
searchDir = fileparts(mfolder); % Get parent folder

% Define the file types to search
filePattern = fullfile(searchDir, '**', '*.m'); % Searches all .m files recursively
files = dir(filePattern);

% Define the keywords to search for
keywords = {k1, k2};

% Store matching files
matchingFiles = {};

% length(files)

for i = 1:length(files)
    filePath = fullfile(files(i).folder, files(i).name);
    
    % Read file content
    try
        fileContent = fileread(filePath);        
        % Check if all keywords exist in the file
        if all(cellfun(@(kw) contains(fileContent, kw), keywords))
            matchingFiles{end+1} = filePath; %#ok<AGROW>
        end
    catch
        warning('Could not read file: %s', filePath);
    end
end

matchingFiles = setdiff(matchingFiles, {[mfilename('fullpath') '.m']});

% Display results
if isempty(matchingFiles)
    disp('No files found containing all keywords.');
else
    disp('Files containing all keywords:');
    disp(matchingFiles);
    if strcmp('Yes', questdlg('Open them?'))
        for k = 1:length(matchingFiles)
            edit(matchingFiles{k})
        end
    end
end


