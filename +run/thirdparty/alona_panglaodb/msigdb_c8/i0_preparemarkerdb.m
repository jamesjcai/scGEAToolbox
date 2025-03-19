species = 'hs';
switch lower(species)
    case {'human', 'hs'}
        url = 'https://scgeatool.github.io/data/msigdb/c8.all.v2022.1.Hs.json';
    case {'mouse', 'mm'}
        url = 'https://scgeatool.github.io/data/msigdb/m8.all.v2022.1.Mm.json';
end

val = webread(url);
setnames = fields(val);

fid = fopen(sprintf('markerlist_%s.txt', species), 'w');

S = [];
for k = 1:length(setnames)
    a = string(val.(setnames{k}).geneSymbols);
    fprintf(fid, '%s\t%s\n', setnames{k}, sprintf('%s,', a));
    S = [S; a]; %#ok
end
fclose(fid);

N = length(S);
t = tabulate(S);
f = cell2mat(t(:, 3));
w = 1 + sqrt((max(f) - f)/(max(f) - min(f)));
genelist = string(t(:, 1));

fid = fopen(sprintf('markerweight_%s.txt', species), 'w');
for k = 1:length(genelist)
    fprintf(fid, '%s\t', genelist(k));
    fprintf(fid, '%f\n', w(k));
end
fclose(fid);
Tm = readtable('markerlist_hs.txt', 'ReadVariableNames', false);
Tw = readtable('markerweight_hs.txt', 'ReadVariableNames', false);
save marker_hs Tm Tw
