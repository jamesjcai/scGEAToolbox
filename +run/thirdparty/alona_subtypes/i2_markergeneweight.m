species='hs';
ctype='neurons';
ctype='tcells';
T=readtable(sprintf('markerlist_%s_%s.txt',species,ctype),'ReadVariableNames',false,'Delimiter','\t');

customf=sprintf('markerlist_%s_%s_custom.txt',species,ctype);
if exist(customf,'file')
    T2=readtable(customf,'ReadVariableNames',false,'Delimiter','\t');
    T=[T;T2];
end

% T.Var1=i_makeuniquename(T.Var1);
% writetable(T,sprintf('markerlist_%s.txt',species),...
%     'WriteVariableNames',false,'Delimiter','\t');

%%
s=upper(string(T.Var2));
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
w=1+sqrt((max(f)-f)/(max(f)-min(f)));
genelist=string(t(:,1));

fid=fopen(sprintf('markerweight_%s_%s.txt',species,ctype),'w');
for k=1:length(genelist)
    fprintf(fid,'%s\t',genelist(k));
    fprintf(fid,'%f\n', w(k));    
end
fclose(fid);

