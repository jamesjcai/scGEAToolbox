species = 'mm';
T = readtable(sprintf('markerlist_%s.txt', species), 'ReadVariableNames', false, 'Delimiter', '\t');
s = string(T.Var2);
S = [];
for k = 1:length(s)
    a = strsplit(s(k), ',');
    a = a(1:end-1);
    S = [S, a];
end
N = length(S);
t = tabulate(S);
f = cell2mat(t(:, 3));
frange = max(f) - min(f);
if frange == 0
    w = ones(size(f));
else
    w = 1 + sqrt((max(f) - f) / frange);
end
genelist = string(t(:, 1));

fid = fopen(sprintf('markerweight_%s.txt', species), 'w');
for k = 1:length(genelist)
    fprintf(fid, '%s\t', genelist(k));
    fprintf(fid, '%f\n', w(k));
end
fclose(fid);