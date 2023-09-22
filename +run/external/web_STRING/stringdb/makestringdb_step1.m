t = readtable('mart_export_mouse.txt');
t(t.ProteinStableID == "", :) = [];
mapObj = containers.Map(t.ProteinStableID, t.GeneName);

T = readtable('10090.protein.links.v11.0.txt');

%%
tic
T.protein1 = extractAfter(T.protein1, 6);
T.protein2 = extractAfter(T.protein2, 6);
i = ismember(T.protein1, t.ProteinStableID);
j = ismember(T.protein2, t.ProteinStableID);
T = T(i & j, :);
T.protein1 = values(mapObj, T.protein1);
T.protein2 = values(mapObj, T.protein2);
toc

writetable(T, 'stringdb_mouse.txt');
