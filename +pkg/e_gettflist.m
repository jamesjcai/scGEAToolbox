function [T] = e_gettflist(speciesid)
% Calcute predefined cell scores (marker list in cellscores.txt/xlsx)
%
% see also: SC_CELLSCORE_UCELL, SC_CELLSCORE_ADMDL, SC_CELLCYCLESCORING

if nargin<1, speciesid = 'hs'; end
pw1 = fileparts(mfilename('fullpath'));
switch lower(speciesid)
    case {'hs', 'human'}
        fname = fullfile(pw1, '..', 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_hs.mat');
    case {'mm', 'mouse'}
        fname = fullfile(pw1, '..', 'assets', 'DoRothEA_TF_Target_DB', 'dorothea_mm.mat');
    otherwise
        error('TF database is not available for the species.');
end
load(fname, 'T');
T = T(T.mor > 0, :);
end
