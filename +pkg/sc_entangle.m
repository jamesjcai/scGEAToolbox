

x=sce.X(ismember(sce.g,{'Mrpl15','Lypla1','Atp6v1h'}),:);
y=(x>0)';
[a,~,b]=unique(y+0,'rows');
n=grpstats(b,b,@numel);
figure;
bar(n)