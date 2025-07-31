T1 = readtable('outfile.txt', 'ReadVariableNames', false, 'Delimiter', '\t');
T1.Var1 = matlab.lang.makeUniqueStrings(T1.Var1);
writetable(T1, 'markerlist.txt', ...
    'WriteVariableNames', false, 'Delimiter', '\t');

%%
s = (string(T1.Var2));
S = [];
for k = 1:length(s)
    a = strsplit(s(k), ',');
    a = strtrim(a);
    if strlength(a(end)) == 0
        a = a(1:end-1);
    end
    S = [S, a];
end

%%
N = length(S);
t = tabulate(S);
f = cell2mat(t(:, 3));
w = 1 + sqrt((max(f) - f)/(max(f) - min(f)));
genelist = string(t(:, 1));

fid = fopen('markerweight.txt', 'w');
for k = 1:length(genelist)
    fprintf(fid, '%s\t', genelist(k));
    fprintf(fid, '%f\n', w(k));
end
fclose(fid);

%%
Tm = readtable('markerlist.txt', 'ReadVariableNames', false);
Tw = readtable('markerweight.txt', 'ReadVariableNames', false);
save marker_hsHPA Tm Tw
delete('markerlist.txt')
delete('markerweight.txt')
delete('outfile.txt')
