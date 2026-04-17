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

% Use native MATLAB umap (Statistics and Machine Learning Toolbox, R2026a+)
% Parameters chosen to match ml_UMAP defaults:
%   NumNeighbors=15  (n_neighbors=15), Distance='euclidean' (metric),
%   EmbeddingDensity=1 (corresponds to min_dist=0.3 in ml_UMAP),
%   NumEpochs=200 (ml_UMAP uses 200 for large / 500 for small datasets),
%   LearnRate=1, Reproducible='on' (randomize=false in ml_UMAP).
%   Initialization='pca' is used; ml_UMAP defaults to 'spectral' but the
%   native function does not support spectral initialization.
if ~isMATLABReleaseOlderThan('R2026a')
    s = umap(full(X).', ...
        NumDimensions=ndim, ...
        NumNeighbors=15, ...
        Distance='euclidean', ...
        EmbeddingDensity=1, ...
        NumEpochs=200, ...
        LearnRate=1, ...
        Reproducible='on');
    return;
end

% Fall back to ml_UMAP for older MATLAB versions
done = false;
try
    s = run.ml_UMAP(X, ndim);
    done = true;
catch ME
    warning(ME.message);
end
if ~done
    disp('Try older version Matlab UMAP (v4.4).');
    s = run.ml_UMAP(X, ndim, 44);
end

end
