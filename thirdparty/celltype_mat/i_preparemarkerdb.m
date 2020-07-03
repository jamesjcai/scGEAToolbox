T=readtable('PanglaoDB_markers_27_Mar_2020.tsv','filetype','text');

a=string(unique(T.cellType));
gt=string(T.cellType);
to=string(T.officialGeneSymbol);
fid=fopen('markerlist_panglaodb.txt','w');
for k=1:length(a)
    k
    idx=find(gt==a(k));
    fprintf(fid,'%s\t',a(k));
    fprintf(fid,'%s,', to(idx));
    fprintf(fid,'\n');
end
fclose(fid);
