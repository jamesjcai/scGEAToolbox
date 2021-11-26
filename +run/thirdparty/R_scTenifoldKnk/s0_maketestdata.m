NCELLS=2000;
NGENES=400;

rng(123);
X=nbinrnd(20,0.98,NGENES,NCELLS);
sce=SingleCellExperiment(X);
if exist('./input.h5','file'), delete('./input.h5'); end

h5create('input.h5', '/X', size(sce.X));
h5write('input.h5', '/X', sce.X);
h5create('input.h5', '/g', size(sce.g),'Datatype','string');
h5write('input.h5', '/g', sce.g);
targetg=sce.g(1);
h5create('input.h5', '/targetg', size(targetg),'Datatype','string');
h5write('input.h5', '/targetg', targetg);
