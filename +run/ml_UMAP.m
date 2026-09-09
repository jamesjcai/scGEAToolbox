function [s] = ml_UMAP(X, ndim, nneighbors)
% ML_UMAP - UMAP embedding using the bundled external/ml_umap45 package
%
%   S = ML_UMAP(X, NDIM) embeds the cells of the gene-by-cell matrix X into
%   NDIM dimensions. NNEIGHBORS defaults to 15, matching the UMAP class.
%
%   This is the legacy fallback used only on MATLAB releases older than
%   R2026a; SC_UMAP calls the built-in UMAP on R2026a and newer.

if nargin < 3, nneighbors = 15; end
if nargin < 2, ndim = 3; end

pw1 = fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)
    % UMAP.m needs both the package root and its util/ subfolder (Args,
    % MatBasics, PopUp, String, ...). Those class names are generic enough
    % to shadow other code, so the entries come off the path on the way out.
    umappth = fullfile(fileparts(pw1), 'external', 'ml_umap45');
    umapcleanup = pkg.i_addpathtemp(umappth, ...
        fullfile(umappth, 'util'));   %#ok<NASGU>
end

data = transpose(X);

ncells = size(data, 1);
if ncells > 500
    if ~(ismcc || isdeployed)
        % svdpca lives in external/ml_PHATE.
        phatecleanup = pkg.i_addpathtemp( ...
            fullfile(fileparts(pw1), 'external', 'ml_PHATE'));   %#ok<NASGU>
    end
    data = svdpca(data, 50, 'random');
end

u = UMAP;
u.n_components = ndim;
u.n_neighbors = nneighbors;
u.setMethod(pkg.i_umapmethod());
s = u.fit_transform(data);
end
