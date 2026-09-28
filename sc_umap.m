function s = sc_umap(X, ndim, donorm, dolog1p)
% UMAP embedding of cells
% s=sc_umap(X,3);

% see also: SC_TSNE, SC_PHATE

if nargin < 4, dolog1p = true; end
if nargin < 3, donorm = true; end
if nargin < 2, ndim = 3; end

if donorm
    X = sc_norm(X, 'type', 'libsize');
    disp('Library-size normalization...done.')
end
if dolog1p
    X = log1p(X);
    disp('log1p transformation...done.')
end

if isMATLABReleaseOlderThan('R2026a')
    % Fall back to the File Exchange UMAP add-on (#71902) on older MATLAB.
    % Note: that package ships MEX binaries for Windows and Intel Mac only.
    % ML_UMAP takes the gene-by-cell matrix and does its own PCA, so hand it
    % the untransposed, unreduced input.
    s = run.ml_UMAP(X, ndim);
    return;
end

data = full(X).';

% Reduce to 50 PCs before the neighbor search, matching RUN.ML_UMAP and the
% NumPCAComponents=50 that SC_TSNE passes to TSNE. The native UMAP has no
% input-reduction option of its own -- Initialization='pca' only seeds the
% starting embedding coordinates -- so the reduction has to happen here.
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

% Use native MATLAB umap (Statistics and Machine Learning Toolbox, R2026a+)
% Parameters chosen to match ml_UMAP defaults:
%   NumNeighbors=15  (n_neighbors=15), Distance='euclidean' (metric),
%   EmbeddingDensity=1 (corresponds to min_dist=0.3 in ml_UMAP),
%   NumEpochs=200 (ml_UMAP uses 200 for large / 500 for small datasets),
%   LearnRate=1, Reproducible='on' (randomize=false in ml_UMAP).
%   Initialization defaults to 'pca'; ml_UMAP defaults to 'spectral' but the
%   native function does not support spectral initialization.
s = umap(data, ...
    NumDimensions=ndim, ...
    NumNeighbors=15, ...
    Distance='euclidean', ...
    EmbeddingDensity=1, ...
    NumEpochs=200, ...
    LearnRate=1, ...
    Reproducible='on');

end
