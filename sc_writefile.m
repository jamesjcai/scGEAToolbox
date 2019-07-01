function sc_writefile(filename,X,genelist)
t=table();
t.genes=string(genelist);
t=[t,array2table(X)];
%[file,path] = uiputfile('data_1.txt');
%if file~=0
writetable(t,filename,'Delimiter','\t','filetype','text');
%end
