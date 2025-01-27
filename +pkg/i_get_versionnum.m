function [v1] = i_get_versionnum
    v1 = '';
    mfolder = fileparts(mfilename('fullpath'));
    vfile = fullfile(mfolder, '..', 'VERSION.mat');

    if exist(vfile, "file")
        data = load(vfile, "v1");
        if isfield(data, "v1") && ~isempty(data.v1)
            v1 = data.v1;
            return;
        end
    end

    tag_version = 'param.version';    
    try
        xfilelocal = fullfile(mfolder,'..', 'scGEAToolbox.prj');
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
    if isempty(v1)
    
        xfilelocal = fullfile(mfolder,'..', 'info.xml');
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
    end
    % Save to VERSION.mat if v1 is non-empty
    if ~isempty(v1) && ~exist(vfile,"file")
        save(vfile, 'v1'); 
    end    
end