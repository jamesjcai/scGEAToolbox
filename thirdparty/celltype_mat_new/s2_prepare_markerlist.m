% T=readtable('clustermole_markers_hs.csv');
T=readtable('clustermole_markers_mm.csv');

C=unique(T.celltype);

%%
fid=fopen('markerlist_mm.txt','w');
% fid=fopen('markerlist_hs.txt','w');
for k=1:length(C)
    k
    idx=string(T.celltype)==C{k};
    fprintf(fid,'%s\t%s',C{k},...
        sprintf('%s,',unique(string(T.gene(idx)))));
    fprintf(fid,'\n');
end
fclose(fid);
