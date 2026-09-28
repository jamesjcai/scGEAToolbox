function [s] = ml_UMAP(X, ndim, nneighbors)
% ML_UMAP - UMAP embedding using the File Exchange UMAP add-on
%
%   S = ML_UMAP(X, NDIM) embeds the cells of the gene-by-cell matrix X into
%   NDIM dimensions. NNEIGHBORS defaults to 15, matching the UMAP class.
%
%   This is the legacy fallback used only on MATLAB releases older than
%   R2026a; SC_UMAP calls the built-in UMAP on R2026a and newer. It needs
%   the "Uniform Manifold Approximation and Projection (UMAP)" add-on
%   (File Exchange #71902); PKG.I_CHECKUMAPADDON says how to install it.

if nargin < 3, nneighbors = 15; end
if nargin < 2, ndim = 3; end

pkg.i_checkumapaddon();

data = transpose(X);

% This is what external/ml_PHATE/svdpca does with method='random', over
% PKG.E_RANDPCA instead of that folder's own RANDPCA. The two give bit-
% identical scores, but e_randPCA puts the caller's random stream back
% afterwards, where RANDPCA leaves the session parked on rng('default').
ncells = size(data, 1);
if ncells > 500
    npc = min(50, size(data, 2));   % RANDPCA refuses k > the smallest dimension
    data = data - mean(data, 1);
    [coeff, ~, ~] = pkg.e_randPCA(data.', npc);
    data = data*coeff;
end

u = UMAP;
u.n_components = ndim;
u.n_neighbors = nneighbors;
u.setMethod(pkg.i_umapmethod());
s = u.fit_transform(data);
end
