function [output] = e_ensembl2symbol(input, species)

if nargin < 2, species = 'human'; end
pw1 = fileparts(mfilename('fullpath'));
%pth = fullfile(pw1, '..', 'assets', 'Biomart', ...
%    sprintf('mart_export_%s.txt', species));
%T = readtable(pth);
pth = fullfile(pw1, '..', 'assets', 'Ensembl', ...
    sprintf('Ensembl_%s_genes.mat', species));

load(pth,"T");

keySet = T.GeneStableID;
valueSet = T.GeneName;
output = input;

idxv = strlength(string(valueSet)) == 0;
valueSet(idxv) = keySet(idxv);

if isstring(input)
    M = containers.Map(string(keySet), string(valueSet));
    %input=string({'ENSG00000121410', 'ENSG00000121411'});
    ix = isKey(M, cellstr(input));
    output(ix) = values(M, cellstr(input(ix)));
elseif iscell(input)
    M = containers.Map(keySet, valueSet);
    %input=string({'ENSG00000121410', 'ENSG00000121411'});
    ix = isKey(M, input);
    output(ix) = values(M, input(ix));
end
