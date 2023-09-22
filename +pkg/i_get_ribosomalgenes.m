function [g] = i_get_ribosomalgenes
%Get ribosomal genes

pw1 = fileparts(mfilename('fullpath'));
ribosomalfile = fullfile(pw1, 'ribosomal.txt');
if ~exist(ribosomalfile, 'file')
    options = weboptions('Timeout', 21);
    % fname=tempname;
    disp('Downloading ribosomal gene names...');
    websave(ribosomalfile, 'https://www.genenames.org/cgi-bin/genegroup/download?id=1054&type=branch', options);
    % t=readtable('a.txt','PreserveVariableNames',true);
end
warning off
t = readtable(ribosomalfile, 'Range', 'A:B');
warning on
g = string(t.ApprovedSymbol);
% delete(fname);
end
