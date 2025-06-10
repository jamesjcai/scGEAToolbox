function [g] = i_get_ribosomalgenes

    pw1 = fileparts(mfilename('fullpath'));
    ribosomalfile = fullfile(pw1, '..','assets', 'HGNC', 'ribosomal.txt');
    if ~exist(ribosomalfile, 'file')    
        options = weboptions('Timeout', 21);
        % fname=tempname;
        disp('Downloading ribosomal gene names...');
        websave(ribosomalfile, 'https://www.genenames.org/cgi-bin/genegroup/download?id=1054&type=branch', options);
        % t=readtable('a.txt','PreserveVariableNames',true);
    end
    warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
    t = readtable(ribosomalfile, 'Range', 'A:B', ...
        'VariableNamingRule', 'modify');
    g = string(t.ApprovedSymbol);
    % delete(fname);
end
