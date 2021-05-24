species='mm';
T1=readtable(sprintf('markerlist_%s_panglaodb.txt',species),'ReadVariableNames',false,'Delimiter','\t');
T2=readtable(sprintf('markerlist_%s_custom.txt',species),'ReadVariableNames',false,'Delimiter','\t');
T=[T1;T2];
s=upper(string(T.Var2));
S=[];
for k=1:length(s)
    a=strsplit(s(k),',');
    a=a(1:end-1);
    S=[S,a];
end
N=length(S);
t=tabulate(S);
f=cell2mat(t(:,3));
w=1+sqrt((max(f)-f)/(max(f)-min(f)));
genelist=string(t(:,1));

fid=fopen(sprintf('markerweight_%s.txt',species),'w');
for k=1:length(genelist)
    fprintf(fid,'%s\t',genelist(k));
    fprintf(fid,'%f\n', w(k));    
end
fclose(fid);

