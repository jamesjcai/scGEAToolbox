function [versionStr] = i_get_versionnum

    versionStr = '';
    mfolder = fileparts(mfilename('fullpath'));
    %{
    tag_version = 'param.version';
    xfilelocal = fullfile(mfolder,'..', 'scGEAToolbox.prj');
    if exist(xfilelocal, 'file')
        try            
            fid = fopen(xfilelocal, 'r');
            C = textscan(fid, '%s', 'delimiter', '\n');
            fclose(fid);
            a = C{1};
            x = a(contains(a, sprintf('<%s>',tag_version)));
            a1 = strfind(x, sprintf('<%s>',tag_version));
            a2 = strfind(x, sprintf('</%s>',tag_version));
            v1 = extractBetween(x, a1{1}+length(sprintf('<%s>',tag_version)), a2{1}-1);
            % v1 = strrep(v1{1}, 'scGEAToolbox ', '');
            v1 = v1{1};
        catch ME
            warning(ME.identifier, 'Error reading project file: %s', ME.message);
        end
    end
    %}
    
    %{
    xfilelocal = fullfile(mfolder, '..', 'info.xml');
    
    fid = fopen(xfilelocal, 'r');
    
    if fid == -1
        warning('Could not open file: %s', xfilelocal);
        return;
    end

    try
        C = textscan(fid, '%s', 'delimiter', '\n');
    catch ME
        warning(ME.identifier, 'Error reading file: %s', ME.message);
    end
    fclose(fid);

    try
        a = C{1};
        % https://www.mathworks.com/matlabcentral/answers/359034-how-do-i-replace-textread-with-textscan
        x = a(contains(a, '<name>'));
        a1 = strfind(x, '<name>');
        a2 = strfind(x, '</name>');
        v1 = extractBetween(x, a1{1}+length('<name>'), a2{1}-1);
        v1 = strrep(v1{1},'scGEAToolbox ','');
    catch ME
        warning(ME.identifier, 'Error reading file content: %s', ME.message);
    end
    %}

    docFolder   = fullfile(mfolder, "..", "doc");
    xfilelocal  = fullfile(docFolder, "helptoc.xml");

    % Resolve to absolute path
    xfilelocal = char(java.io.File(xfilelocal).getCanonicalPath());

    % xfilelocal = fullfile(mfolder, '..', 'doc', 'helptoc.xml');
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