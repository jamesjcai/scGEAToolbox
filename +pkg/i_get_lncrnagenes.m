function [g] = i_get_lncrnagenes

    pw1 = fileparts(mfilename('fullpath'));
    lncrnafile = fullfile(pw1, '..','resources', 'HGNC', 'lncrna.txt');
    if ~exist(lncrnafile, 'file')
        % options = weboptions('Timeout', 21);
        % % fname=tempname;
        % disp('Downloading ribosomal gene names...');
        % websave(lncrnafile, 'https://www.genenames.org/cgi-bin/genegroup/download?id=1054&type=branch', options);
        % % t=readtable('a.txt','PreserveVariableNames',true);
    end
    warning off
    t = readtable(lncrnafile,'ReadVariableNames',false);
    warning on
    g = string(t.Var1);
    % delete(fname);
end
