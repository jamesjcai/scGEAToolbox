species = 'hs';
T1 = readtable(sprintf('markerlist_%s_panglaodb.txt', species), 'ReadVariableNames', false, 'Delimiter', '\t');
T2 = readtable(sprintf('markerlist_%s_custom.txt', species), 'ReadVariableNames', false, 'Delimiter', '\t');
% T3 and T4 downloaded from Enrichr
% T3=readtable('Descartes_Cell_Types_and_Tissue_2021.txt','ReadVariableNames',false,'Delimiter','\t');
%T4=readtable('CellMarker_Augmented_2021.txt','ReadVariableNames',false,'Delimiter','\t');

% tmpT=readtable('CellMarkerAugmented2021_ori.txt');
% switch species
%     case 'hs'
%         tmpT=tmpT(tmpT.speciesType=="Human",:);
%     case 'mm'
%         tmpT=tmpT(tmpT.speciesType=="Mouse",:);
% end
% T5=table(tmpT.cellName,tmpT.geneSymbol);

switch species
    case 'hs'
        T = [T1; T2];
    case 'mm'
        T = [T1; T2];
end
T.Var1 = i_makeuniquename(T.Var1);
writetable(T, sprintf('markerlist_%s.txt', species), ...
    'WriteVariableNames', false, 'Delimiter', '\t');

T(T.Var1 == "Gamma delta T cells", :) = [];

%%
s = upper(string(T.Var2));
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

fid = fopen(sprintf('markerweight_%s.txt', species), 'w');
for k = 1:length(genelist)
    fprintf(fid, '%s\t', genelist(k));
    fprintf(fid, '%f\n', w(k));
end
fclose(fid);


Tm = T;
Tw = readtable(sprintf('markerweight_%s.txt', species));
save(sprintf('marker_%s', species), 'Tw', 'Tm');
