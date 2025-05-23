function [matchingFiles] = e_openallfiles(k1, k2, bequiet)

if nargin < 3, bequiet = false; end
% if nargin<1, k1='parentfig'; k2='myQuestdlg'; end
% if nargin<1, k1='FigureHandle'; k2='inputdlg'; end
%if nargin<1, k1='isempty(answer'; k2='myQuestdlg'; end
%if nargin<1, k1='isempty'; k2='switch'; end
%if nargin<1, k1='InitialValue'; k2='myListdlg'; end
if nargin<1, k1='readtable'; k2='warn'; end



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
        % if strcmp(filePath, 'D:\GitHub\scGEAToolbox\+pkg\i_makestandalone.m')
        %     disp('ok')
        %     fileContent = fileread(filePath);
        %     assignin("base","filecontent", fileContent);
        % end
        fileContent = fileread(filePath);
        % Check if all keywords exist in the file
        if all(cellfun(@(kw) contains(fileContent, kw), keywords))
        %if ~contains(fileContent, keywords{1}) && contains(fileContent, keywords{2})
            matchingFiles{end+1} = filePath; %#ok<AGROW>
        end
    catch
        warning('Could not read file: %s', filePath);
    end
end

matchingFiles = setdiff(matchingFiles, {[mfilename('fullpath') '.m']});

% Display results
if isempty(matchingFiles)
    % disp('No files found containing all keywords.');
else
    %disp('Files containing all keywords:');
    %disp(matchingFiles);
    if ~bequiet
    if strcmp('Yes', questdlg('Open them?'))
        for k = 1:length(matchingFiles)
            edit(matchingFiles{k})
        end
    end
    end
end


