function [v1] = i_getversionnum

    mfolder = fileparts(mfilename('fullpath'));
    xfile = 'scGEAToolbox.prj';
    tag = 'param.version';
    
    v1 = '';
    try
        xfilelocal = fullfile(mfolder,'..', xfile);
        fid = fopen(xfilelocal, 'r');
        C = textscan(fid, '%s', 'delimiter', '\n');
        fclose(fid);
        a = C{1};
        x = a(contains(a, sprintf('<%s>',tag)));
        a1 = strfind(x, sprintf('<%s>',tag));
        a2 = strfind(x, sprintf('</%s>',tag));
        v1 = extractBetween(x, a1{1}+length(sprintf('<%s>',tag)), a2{1}-1);
        % v1 = strrep(v1{1}, 'scGEAToolbox ', '');
        v1 = v1{1};
    catch
    end
    if isempty(v1)
    try
        xfilelocal = fullfile(mfolder,'..', 'info.xml');
        fid = fopen(xfilelocal, 'r');
        C = textscan(fid, '%s', 'delimiter', '\n');
        fclose(fid);
        a = C{1};
        % https://www.mathworks.com/matlabcentral/answers/359034-how-do-i-replace-textread-with-textscan
        x = a(contains(a, '<name>'));
        a1 = strfind(x, '<name>');
        a2 = strfind(x, '</name>');
        v1 = extractBetween(x, a1{1}+length('<name>'), a2{1}-1);
        v1 = strrep(v1{1},'scGEAToolbox ','');
    catch
    end
    end
end