function [sce] = r_seurat(X, genelist, wkdir, isdebug)

if nargin < 3, wkdir = tempdir; end
if nargin < 4, isdebug = false; end

oldpth = pwd();
[isok, msg, codepath] = commoncheck_R('R_Seurat');
if ~isok, error(msg); end

if ~isok, error(msg);
    sce = [];
    return;
end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end



if isa(X, 'SingleCellExperiment') && isnumeric(genelist)
    sce = X;
    ndim = genelist;
else
    if nargin < 2, error("[sce]=run.SeuratWorkflow(X,genelist)"); end
    sce = SingleCellExperiment(X, genelist);
    ndim = 2;
end

tmpfilelist = {'input.mat', 'output.h5', 'g.txt'};
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end

X = sce.X;
if issparse(X), X = full(X); end
save('input.mat', 'X', 'ndim', '-v7.3');
writematrix(sce.g, 'g.txt');

Rpath = getpref('scgeatoolbox', 'rexecutablepath');

codefullpath = fullfile(codepath,'script.R');
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.i_runrcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    s_tsne = h5read('output.h5', '/s_tsne');
    s_umap = h5read('output.h5', '/s_umap');
    c_ident = h5read('output.h5', '/c_ident');
    %X = h5read('output.h5', '/X');
    [c, ~] = findgroups(c_ident);
    sce.c_cluster_id = c;
    sce.struct_cell_clusterings.seurat = c_ident;
    sce.c = c;
    %sce.X = X;

    if ~isfield(sce.struct_cell_embeddings,'umap3d')
        sce.struct_cell_embeddings = setfield(sce.struct_cell_embeddings, 'umap3d', []);
    end
    if ~isfield(sce.struct_cell_embeddings,'umap2d')
        sce.struct_cell_embeddings = setfield(sce.struct_cell_embeddings, 'umap2d', []);
    end
    if ~isfield(sce.struct_cell_embeddings,'tsne3d')
        sce.struct_cell_embeddings = setfield(sce.struct_cell_embeddings, 'tsne3d', []);
    end
    if ~isfield(sce.struct_cell_embeddings,'tsne2d')
        sce.struct_cell_embeddings = setfield(sce.struct_cell_embeddings, 'tsne2d', []);
    end
    

    if size(s_umap,2) == 3
        sce.struct_cell_embeddings.umap3d = s_umap;
    else
        sce.struct_cell_embeddings.umap2d = s_umap;
    end

    if size(s_tsne,2) == 3
        sce.struct_cell_embeddings.tsne3d = s_tsne;
    else
        sce.struct_cell_embeddings.tsne2d = s_tsne;
    end

    sce.s = s_tsne;
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
cd(oldpth);
end
