function [Y, S] = ml_metaviz(X, ndim, showwaitbar, dophate)

if nargin < 4 || isempty(dophate), dophate = true; end
if nargin < 3 || isempty(showwaitbar), showwaitbar = true; end
if nargin < 2 || isempty(ndim), ndim = 2; end
S = [];
pw1 = fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)
    % phate lives in external/ml_PHATE.
    phatecleanup = pkg.i_addpathtemp( ...
        fullfile(fileparts(pw1), 'external', 'ml_PHATE'));   %#ok<NASGU>
end
if isMATLABReleaseOlderThan('R2026a')
    % Fail now, not after the t-SNE and PHATE runs, if the UMAP add-on the
    % pre-R2026a path needs is missing.
    pkg.i_checkumapaddon();
end

nstep = 6 + 1;
usingmmfile = false;
try
    zeros(size(X, 2), size(X, 2), 14, 'single');
catch ME
    disp(ME.message);
    usingmmfile = true;
    disp('Using memory mapping file.');
end

if showwaitbar, fw = gui.myWaitbar([]); end
Xn = log1p(sc_norm(X))';

% Up to 300 components, dropping to 50 when the data is too small to
% support that -- which is what the try/catch here used to accomplish, by
% letting RANDPCA throw its "k must be <= the smallest dimension" error and
% catching it. Asking first is cheaper and does not use an exception for
% flow control. The final MIN only bites where the old catch branch threw
% in its turn and took the whole function down with it.
npc = 300;
if npc > min(size(Xn)), npc = 50; end
npc = min(npc, min(size(Xn)));

% This is what external/ml_PHATE/svdpca did with method='random', over
% PKG.E_RANDPCA instead of that folder's own RANDPCA. The two give bit-
% identical scores, but e_randPCA puts the caller's random stream back
% afterwards, where RANDPCA leaves the session parked on rng('default') --
% here that reset landed immediately before this function's own t-SNE,
% UMAP and PHATE views.
data = Xn - mean(Xn, 1);
[coeff, ~, ~] = pkg.e_randPCA(data.', npc);
data = data*coeff;

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - PCA...', 1/nstep); end
[~, S{1}] = pca(data, NumComponents = ndim);

try
    d = dot(data, data, 2);
    DS = d + d' - 2 * (data * data');
    % [1] Albanie, Samuel. Euclidean Distance Matrix Trick. June, 2019. Available at https://www.robots.ox.ac.uk/%7Ealbanie/notes/Euclidean_distance_trick.pdf.
catch
    DS = pdist2(data, data).^2;
end

try
    if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - MDS...', 1/nstep); end
    S{end+1} = pkg.e_embedbyd(sqrt(DS), ndim, 2);
catch
    % MDS is optional; the other embeddings (KPCA, t-SNE, ...) still run
end

% The SECOND output. PKG.E_KPCA's first output is COEFF, the projection
% coefficients for mapping new data, not an embedding: EIGS returns
% unit-norm eigenvectors and score = Kc*coeff = coeff*latent, so taking
% coeff gives every kernel component equal weight and throws the spectrum
% away. Its own header describes coeff as something you multiply a kernel
% by to obtain an embedding.
if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - KPCA1...', 2/nstep); end
[~, kpcaScore] = pkg.e_kpca(DS, ndim, 30, true);
S{end+1} = kpcaScore;

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - KPCA2...', 2/nstep); end
[~, kpcaScore] = pkg.e_kpca(DS, ndim, 40, true);
S{end+1} = kpcaScore;

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - KPCA3...', 2/nstep); end
[~, kpcaScore] = pkg.e_kpca(DS, ndim, 50, true);
S{end+1} = kpcaScore;

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - TSNE 1/3...', 3/nstep); end
S{end+1} = tsne(data, Perplexity = 30, NumDimensions = ndim);

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - TSNE 2/3...', 3/nstep); end
S{end+1} = tsne(data, Perplexity = 15, NumDimensions = ndim);

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - TSNE 3/3...', 3/nstep); end
S{end+1} = tsne(data, Perplexity = 50, NumDimensions = ndim);

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - UMAP 1/3...', 4/nstep); end
if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - UMAP 2/3...', 4/nstep); end
if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - UMAP 3/3...', 4/nstep); end
if ~isMATLABReleaseOlderThan('R2026a')
    S{end+1} = umap(full(data), NumDimensions=ndim, NumNeighbors=15);
    S{end+1} = umap(full(data), NumDimensions=ndim, NumNeighbors=30);
    S{end+1} = umap(full(data), NumDimensions=ndim, NumNeighbors=50);
else
    S{end+1} = i_legacyumap(data, ndim, 15);
    S{end+1} = i_legacyumap(data, ndim, 30);
    S{end+1} = i_legacyumap(data, ndim, 50);
end

if dophate
    % PHATE finishes with metric MDS, and MDSCALE gives up on some inputs
    % with "Unable to decrease criterion along line search direction". That
    % killed the whole function, including the dozen embeddings already
    % computed, and DOPHATE is on by default -- on three of three simulated
    % datasets RUN.ML_METAVIZ(X) failed outright while RUN.ML_PHATE on the
    % same counts was fine. A meta-visualisation is a consensus over many
    % views, so a view that will not converge should drop out of the
    % consensus rather than end it.
    phateArgs = {{'t', 20, 'ndim', ndim, 'k', 5}, ...
        {'t', 20, 'ndim', ndim, 'k', 15, 'pot_method', 'sqrt'}, ...
        {'t', 20, 'ndim', ndim, 'k', 30}};
    for phateStep = 1:numel(phateArgs)
        if showwaitbar
            gui.myWaitbar([], fw, [], [], sprintf( ...
                'Meta Visualization - PHATE %d/3...', phateStep), 5/nstep);
        end
        try
            S{end+1} = i_guardedphate(sqrt(Xn), phateArgs{phateStep}); %#ok<AGROW>
        catch ME
            warning('ml_metaviz:phateFailed', ...
                'PHATE view %d of 3 did not converge (%s); it is left out.', ...
                phateStep, ME.message);
        end
    end
end

if showwaitbar, gui.myWaitbar([], fw, [], [], 'Meta Visualization - METAVIZ', 6/nstep); end
if usingmmfile
    [Y] = metaviz_memmap(S, ndim);
else
    [Y] = metaviz_tensor(S, ndim);
end
if showwaitbar, gui.myWaitbar([], fw); end
end

function Y = i_guardedphate(data, args)
% Call PHATE and leave the session state as it was found.
%
% PHATE reaches external/ml_PHATE/randPCA, which opens its randomized
% branch with rng('default') and a bare "warning off" and undoes neither.
% One default RUN.ML_METAVIZ call therefore parked the whole session on
% seed 0 and silenced every later warning. RUN.ML_PHATE guards its own
% PHATE call the same way; this is the other call site. The folder is
% third-party and out of scope for the repository's hygiene rules, so the
% restoring belongs on this side of the call.
%
% A wrapper rather than a save/restore pair around the loop in the caller:
% the ONCLEANUP objects fire the moment this function returns, so the state
% is already back when the caller's CATCH raises ml_metaviz:phateFailed --
% which PHATE's own "warning off" would otherwise have suppressed for every
% view after the first -- and it holds when PHATE throws.
rngState = rng();
restoreRng = onCleanup(@() rng(rngState));
warnState = warning();
restoreWarn = onCleanup(@() warning(warnState));
Y = phate(data, args{:});
end

function reduction = i_legacyumap(data, ndim, nneighbors)
% Embed an already normalized/reduced cells-by-features matrix with the
% File Exchange UMAP add-on (pre-R2026a fallback only).
u = UMAP;
u.n_components = ndim;
u.n_neighbors = nneighbors;
u.verbose = false;
u.setMethod(pkg.i_umapmethod());
reduction = u.fit_transform(data);
end

% PCA: the fast SVD function svds from R package rARPACK with embedding dimension k=2.
% MDS: the basic R function cmdscale with embedding dimension k=2.
% Sammon: the R function sammon from R package MASS with embedding dimension k=2.
% LLE: the R function lle from R package lle with parameters m=2, k=20, reg=2.
% HLLE: the R function embed from R package dimRed with parameters method="HLLE",knn=20, ndim=2.
% Isomap: the R function embed from R package dimRed with parameters method="Isomap",knn=20, ndim=2.
% kPCA1&2: the R function embed from R package dimRed with parameters method="kPCA",kpar=list(sigma=width), ndim=2, where we set width=0.01 for kPCA1 and width=0.001 for kPCA2.
% LEIM: the R function embed from R package dimRed with parameters ndim=2 and method ="LaplacianEigenmaps".
% UMAP1&2: the R function umap from R package uwot with parameters n neighbors=n,n components=2, where we set n=30 for UMAP1 and width=50 for UMAP2.
% tSNE1&2: the R function embed from R package dimRed with parameters method="tSNE",perplexity=n, ndim=2, where we set n=10 for tSNE1 and n=50 for tSNE2.
% PHATE1&2: the R function phate from R package phateR with parameters knn=n, ndim=2,where we set n=30 for PHATE1 and n=50 for PHATE2.


function [Y] = metaviz_memmap(Sinput, ndim, methodid)
if nargin < 3, methodid = 1; end
if nargin < 2, ndim = 2; end

K = length(Sinput); % K = number of embeddings
n = size(Sinput{1}, 1); % n = number of cells
w = zeros(K, n);

%%


mmf = tempname;
fileID = fopen(mmf, 'w');
for k = 1:K
    fwrite(fileID, zeros([n * n, 1]), 'single');
end
fclose(fileID);
m = memmapfile(mmf, 'Format', 'single', 'Writable', true);

for k = 1:K
    d = pdist2(Sinput{k}, Sinput{k});
    m.Data((n^2)*(k - 1)+1:(n^2)*(k)) = d ./ vecnorm(d);
end

%%
m.Offset = 0;
for x = 1:n % n of cells
    % Column x of each n-by-n block, matching the tensor path. This was
    % reshape(m.Data(x:n:N), [n, K]), and striding by n from x walks ROW x
    % of a column-major block -- the same wrong slice, for the same
    % reason. The blocks are laid out column-major, so column x of block k
    % starts at (x-1)*n + 1 + (n*n)*(k-1), which is exactly the indexing
    % the combination loop below already uses.
    d = zeros(n, K, 'single');
    for k = 1:K
        s = (x - 1)*n + 1 + (n*n)*(k - 1);
        d(:, k) = m.Data(s:s + n - 1);
    end
    S = 1 - squareform(pdist(d', 'cosine'));

    [v, ~] = eigs(double(S), 1);
    w(:, x) = abs(v);
end

m.Offset = 0;
M = zeros(n, n);
for l = 1:n % cell
    d = zeros(n, 1);
    for k = 1:K % type of embedding
        s1 = (l - 1) * n + 1;
        s2 = (n * n) * (k - 1);
        s = s1 + s2;
        t = s + n - 1;
        d = d + w(k, l)' .* m.Data(s:t);
    end
    M(:, l) = d;
end
M = 0.5 * (M + M.');

if exist(mmf, 'file') == 2, delete(mmf); end
[Y] = pkg.e_embedbyd(M, ndim, methodid);

end


function [Y] = metaviz_tensor(Sinput, ndim, methodid)

if nargin < 3, methodid = 1; end
if nargin < 2, ndim = 2; end

K = length(Sinput); % K = number of embeddings
n = size(Sinput{1}, 1); % n = number of cells
w = zeros(K, n);

%%
D = zeros(n, n, K, 'single');
for k = 1:K
    d = pdist2(Sinput{k}, Sinput{k});
    D(:, :, k) = d ./ vecnorm(d);
    % (n^2)*(k-1)+1:(n^2)*(k)
end

%%
% D(:, x, :), not D(x, :, :). The line above normalises COLUMNS to unit
% norm, so column x of each slice is cell x's distance profile scaled by
% one constant, while row x has every entry divided by a different
% column's norm. Cosine distance ignores a per-vector constant, so the
% correct slice gives exactly the similarities the raw distances would --
% checked at 3e-16 -- whereas the row slice differs from them by ~7e-4.
% The combination loop below already reads D(:, i, k), so the weights were
% computed on one slicing and applied to another.
for x = 1:n % n of cells
    S = 1 - squareform(pdist(squeeze(D(:, x, :))', 'cosine'));
    [v, ~] = eigs(double(S), 1);
    w(:, x) = abs(v);
end

M = zeros(n, n);
for i = 1:n % cell
    d = zeros(n, 1);
    for k = 1:K % type of embedding
        d = d + w(k, i)' .* D(:, i, k);
    end
    M(:, i) = d;
end
M = 0.5 * (M + M.');

[Y] = pkg.e_embedbyd(M, ndim, methodid);

end
