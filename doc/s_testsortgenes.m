T=readtable('genelist.txt');
c=T.Chromosome_scaffoldName;

[y,idx]=ismember(g,string(T.GeneName));
[~,idx2]=sort(idx(y));


%%
c1=c(y,:);
c1=c1(idx2,:);
c2=c(~y);
c2(:)=23;
cnew=[c1;c2];

X1=X(y,:);
X1=X1(idx2,:);
X2=X(~y,:);
Xnew=[X1;X2];

g1=g(y,:);
g1=g1(idx2);
g2=g(~y);
gnew=[g1;g2];

%%
A=(log10(1+log10(1+log10(1+log10(1+zscore(Xnew,0,'all'))))));
[~,idx3]=sort(sce.c_cluster_id);

imagesc(A(:,idx3))
p=cumsum(find(diff(cnew)>0));
for k=1:length(p)
    yline(p(k),'-w',sprintf('%d',k))
end
