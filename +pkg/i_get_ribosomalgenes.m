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
% Restored on the way out. The bare warning('off', ...) left this
% identifier disabled for the rest of the MATLAB session, and this
% function is on ordinary paths -- ten.sctenifoldnet and
% SingleCellExperiment.rmribosomalgenes both call it -- so merely filtering
% ribosomal genes silenced the warning everywhere afterwards.
% GUI.I_READTABLEFILE already does it this way for the same identifier.
warnstate = warning('off', 'MATLAB:table:ModifiedAndSavedVarnames');
cleanupObj = onCleanup(@() warning(warnstate));

t = readtable(ribosomalfile, 'Range', 'A:B', ...
'VariableNamingRule', 'modify');
g = string(t.ApprovedSymbol);
% delete(fname);
end
