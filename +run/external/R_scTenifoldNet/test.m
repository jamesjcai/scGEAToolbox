NCELLS=2000;
NGENES=400;

rng(123);
X=nbinrnd(20,0.98,NGENES,NCELLS);
sce1=SingleCellExperiment(X);

rng(223);
X=nbinrnd(20,0.98,NGENES,NCELLS);
sce2=SingleCellExperiment(X);

[T]=run.r_sctenifoldnet(sce1,sce2);

