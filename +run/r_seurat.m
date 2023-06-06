function [sce]=r_seurat(X,genelist)

isdebug=false;
oldpth=pwd();
[isok,msg]=commoncheck_R('R_Seurat');
if ~isok, error(msg); end

if ~isok, error(msg); sce=[]; return; end
if isa(X,'SingleCellExperiment') && isnumeric(genelist)
    sce=X;
    ndim=genelist;
else
    if nargin<2, error("[sce]=run.SeuratWorkflow(X,genelist)"); end
    sce=SingleCellExperiment(X,genelist);
    ndim=2;
end

tmpfilelist={'input.mat','output.h5','g.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

X=sce.X;
save('input.mat','X','ndim','-v7.3');
writematrix(sce.g,'g.txt');

Rpath=getpref('scgeatoolbox','rexecutablepath');
pkg.RunRcode('script.R',Rpath);

if exist('./output.h5','file')
    s_tsne=h5read('output.h5','/s_tsne');
    s_umap=h5read('output.h5','/s_umap');
    c_ident=h5read('output.h5','/c_ident');
    [c,~]=grp2idx(c_ident);
    sce.c_cluster_id=c;
    sce.struct_cell_clusterings.seurat=c_ident;
    sce.c=c;
    sce.struct_cell_embeddings.umap=s_umap;
    sce.struct_cell_embeddings.tsne=s_tsne;
    sce.s=s_tsne;    
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
