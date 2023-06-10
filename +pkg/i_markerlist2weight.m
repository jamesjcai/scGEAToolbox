function [Tm,Tw]=i_markerlist2weight(sce)
%Tm=readtable('markerlist_hs.txt','ReadVariableNames',false);
%celltypev=string(Tm.Var1);
%markergenev=string(Tm.Var2);
Tm=[]; Tw=[];

if nargin<1, sce=[]; end
if isempty(sce)
    indata=sprintf('Cell type 1\tgene1,gene2,gene3\nCell type 2\tgene4,gene5');
else
    a=sce.g(randperm(length(sce.g)));
    a1=sprintf('%s,%s,%s',a(1),a(2),a(3));
    a2=sprintf('%s,%s',a(4),a(5));
    indata=sprintf('Cell type 1\t%s\nCell type 2\t%s',a1,a2);
end

% indata=gui.i_getsctypemarkers;

a=inputdlg(sprintf('Format:\nCell type name [TAB] Gene1,Gene2'), ...
    'Markers Input',[10 50],{char(indata)},'on');
if isempty(a), return; end
b=strtrim(string(a{1}));
[c,d]=strtok(b,sprintf('\t'));
d=upper(d);
d=strrep(d,' ','');
Tm=table(strtrim(c),strtrim(d));


s=upper(string(Tm.Var2));
S=[];
for k=1:length(s)
    a=strsplit(s(k),',');
    a=strtrim(a);    
    if strlength(a(end))==0 || isempty(a(end))
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


