function [A, g] = e_prunenet(A, g, whitelist, species)

if nargin < 3, whitelist = []; end
if nargin < 4, species = 'human'; end

%url='https://www.grnpedia.org/trrust/data/trrust_rawdata.mouse.tsv';
url = sprintf('https://www.grnpedia.org/trrust/data/trrust_rawdata.%s.tsv', ...
    lower(species));
T = readtable(url, 'FileType', 'text');
validg = unique([upper(string(T.Var1)); upper(string(T.Var2)); upper(whitelist)]);
y = ismember(upper(g), validg);
A = A(y, y);
g = g(y);
