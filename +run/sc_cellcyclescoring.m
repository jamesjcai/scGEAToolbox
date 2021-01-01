% https://rdrr.io/cran/Seurat/src/R/utilities.R
% https://satijalab.org/seurat/v2.4/cell_cycle_vignette.html

[~,sgenes,g2mgenes]=pkg.i_get_cellcyclegenes;
g=upper(genelist);

[~,idx1]=intersect(g,sgenes);
[~,idx2]=intersect(g,g2mgenes);
i1=randsample(length(g),length(idx1));
i2=randsample(length(g),length(idx2));

%%
score_s=zeros(size(X,2),1);
score_g=zeros(size(X,2),1);
tic
for k=1:size(X,2)
    x=X(:,k);
    score_s(k)=mean(x(idx1))-mean(x(i1));
    score_g(k)=mean(x(idx2))-mean(x(i2));
end
toc