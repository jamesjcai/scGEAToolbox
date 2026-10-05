function [sce] = r_seurat(X, genelist, wkdir, isdebug)

if nargin < 3, wkdir = pkg.i_tempdirfile(); end
if nargin < 4, isdebug = false; end

oldpth = pwd();
cleanupCwd = onCleanup(@() cd(oldpth));
[isok, msg, codepath] = commoncheck_R('R_Seurat');
if ~isok, error('%s', msg); end
if ~isempty(wkdir) && isfolder(wkdir), cd(wkdir); end


givensce = isa(X, 'SingleCellExperiment') && isnumeric(genelist);
if givensce
    sce = X;
    ndim = genelist;
else
    if nargin < 2, error("[sce]=run.r_seurat(X,genelist)"); end
    sce = SingleCellExperiment(X, genelist);
    ndim = 2;
end

tmpfilelist = {'input.mat', 'output.h5', 'g.txt'};
pkg.i_deletefiles(tmpfilelist);   % always clear stale files, so a failed
% run cannot leave a previous run's output to be picked up as this one's

X = sce.X;
if issparse(X), X = full(X); end
save('input.mat', 'X', 'ndim', '-v7.3');
writematrix(sce.g, 'g.txt');

Rpath = getpref('scgeatoolbox', 'rexecutablepath', []);
if isempty(Rpath)
    error('R environment has not been set up.');
end

codefullpath = fullfile(codepath,'script.R');
pkg.i_addwd2script(codefullpath, wkdir, 'R');
pkg.i_runrcode(codefullpath, Rpath);

if exist('output.h5', 'file')
    s_tsne = h5read('output.h5', '/s_tsne');
    s_umap = h5read('output.h5', '/s_umap');
    c_ident = h5read('output.h5', '/c_ident');
    [c, ~] = findgroups(c_ident);
    % Seurat's clusters are stored, not applied, when the caller hands in
    % a dataset: its own clusters stay, and Cluster > Cluster Cells >
    % "Seurat FindClusters" applies these on request. They used to replace
    % C_CLUSTER_ID every time the embedding was run. A dataset made here
    % from a bare matrix has no clusters of its own, so it gets them.
    sce.struct_cell_clusterings.seurat = c;
    if ~givensce
        sce.c_cluster_id = c;
    end
    sce.c = c;

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
else
    error('run.r_seurat:noOutput', ...
        ['R finished but did not write %s to %s. The R console ', ...
         'output above should say why.'], 'output.h5', pwd);
end
if ~isdebug, pkg.i_deletefiles(tmpfilelist); end
end
