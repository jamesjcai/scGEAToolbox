% http://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp

fname = 'c2.cp.v2022.1.Hs.json'; 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);


a=fields(val);
idx=gui.i_selmultidlg(a);

string(val.(a{idx}).geneSymbols)

%%

[indx,tf] = listdlg('PromptString',...
    {'Select one:'},...
     'SelectionMode','single','ListString',a,'ListSize',[300 300]);

