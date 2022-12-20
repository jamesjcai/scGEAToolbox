function [Tm,Tw]=i_markerlist2weight
%Tm=readtable('markerlist_hs.txt','ReadVariableNames',false);
%celltypev=string(Tm.Var1);
%markergenev=string(Tm.Var2);


indata=sprintf('cell type 1\tgene1,gene2,gene3\ncell type 2\tgene4,gene5');
a=inputdlg(sprintf('Format:\nCell type name [TAB] Gene1,Gene2'), ...
    'Markers Input',[10 50],{char(indata)});
if isempty(a), return; end
b=strtrim(string(a{1}));
[c,d]=strtok(b,sprintf('\t'));
Tm=table(strtrim(c),strtrim(d));


s=upper(string(Tm.Var2));
S=[];
for k=1:length(s)
    a=strsplit(s(k),',');
    a=strtrim(a);
    if strlength(a(end))==0
        a=a(1:end-1);
    end
    S=[S,a];
end
%%
N=length(S);
t=tabulate(S);
f=cell2mat(t(:,3));
if max(f)-min(f)<eps
    w=ones(N,1);
else
    w=1+sqrt((max(f)-f)/(max(f)-min(f)));
end
genelist=string(t(:,1));
Tw=table(genelist,w);
Tw.Properties.VariableNames={'Var1','Var2'};


