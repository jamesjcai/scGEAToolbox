function [Y, S] = ml_metaviz(X, ndim, showwaitbar, dophate)

if nargin < 4 || isempty(dophate), dophate = true; end
if nargin < 3 || isempty(showwaitbar), showwaitbar = true; end
if nargin < 2 || isempty(ndim), ndim = 2; end
S = [];
pw1 = fileparts(mfilename('fullpath'));
if ~(ismcc || isdeployed)
    pth = fullfile(pw1, '..', 'external', 'ml_PHATE'); % for calling randmds.m
    addpath(pth);
    pth1 = fullfile(pw1, '..', 'external', 'ml_cbrewer');
    addpath(pth1);
    umapversion = 'ml_umap45';
    pth1 = fullfile(pw1, '..', 'external', umapversion);
    addpath(pth1);
    % pth3 = fullfile(pw1, '..', 'external', 'ml_UMAP', 'umap.jar');
    % javaaddpath(pth3);
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

if showwaitbar, fw = gui.gui_waitbar_adv; end
Xn = log1p(sc_norm(X))';
try
    data = svdpca(Xn, 300, 'random');
catch
    data = svdpca(Xn, 50, 'random');
end
%data=Xn;

if showwaitbar, gui.gui_waitbar_adv(fw, 1/nstep, 'Meta Visualization - PCA...'); end
[~, S{1}] = pca(data, NumComponents = ndim);

try
    d = dot(data, data, 2);
    DS = d + d' - 2 * (data * data');
    % [1] Albanie, Samuel. Euclidean Distance Matrix Trick. June, 2019. Available at https://www.robots.ox.ac.uk/%7Ealbanie/notes/Euclidean_distance_trick.pdf.
catch
    DS = pdist2(data, data).^2;
end

try
    if showwaitbar, gui.gui_waitbar_adv(fw, 1/nstep, 'Meta Visualization - MDS...'); end
    %D=squareform(pdist(data));
    S{end+1} = pkg.e_embedbyd(sqrt(DS), ndim, 2);
catch

end

%[y]=pkg.isomap(log(sc_norm(X)+1)');

if showwaitbar, gui.gui_waitbar_adv(fw, 2/nstep, 'Meta Visualization - KPCA1...'); end
S{end+1} = pkg.kpca(DS, ndim, 30, true);

if showwaitbar, gui.gui_waitbar_adv(fw, 2/nstep, 'Meta Visualization - KPCA2...'); end
S{end+1} = pkg.kpca(DS, ndim, 40, true);

if showwaitbar, gui.gui_waitbar_adv(fw, 2/nstep, 'Meta Visualization - KPCA3...'); end
S{end+1} = pkg.kpca(DS, ndim, 50, true);

if showwaitbar, gui.gui_waitbar_adv(fw, 3/nstep, 'Meta Visualization - TSNE 1/3...'); end
S{end+1} = tsne(data, Perplexity = 30, NumDimensions = ndim);

if showwaitbar, gui.gui_waitbar_adv(fw, 3/nstep, 'Meta Visualization - TSNE 2/3...'); end
S{end+1} = tsne(data, Perplexity = 15, NumDimensions = ndim);

if showwaitbar, gui.gui_waitbar_adv(fw, 3/nstep, 'Meta Visualization - TSNE 3/3...'); end
S{end+1} = tsne(data, Perplexity = 50, NumDimensions = ndim);

if showwaitbar, gui.gui_waitbar_adv(fw, 4/nstep, 'Meta Visualization - UMAP 1/3...'); end
S{end+1} = run_umap_lite(data, 'n_components', ndim, ...
    'n_neighbors', 15, 'verbose', 'none');

if showwaitbar, gui.gui_waitbar_adv(fw, 4/nstep, 'Meta Visualization - UMAP 2/3...'); end
S{end+1} = run_umap_lite(data, 'n_components', ndim, ...
    'n_neighbors', 30, 'verbose', 'none');

if showwaitbar, gui.gui_waitbar_adv(fw, 4/nstep, 'Meta Visualization - UMAP 3/3...'); end
S{end+1} = run_umap_lite(data, 'n_components', ndim, ...
    'n_neighbors', 50, 'verbose', 'none');

if dophate
    if showwaitbar, gui.gui_waitbar_adv(fw, 5/nstep, 'Meta Visualization - PHATE 1/3...'); end
    S{end+1} = phate(sqrt(Xn), 't', 20, 'ndim', ndim, 'k', 5);
    
    if showwaitbar, gui.gui_waitbar_adv(fw, 5/nstep, 'Meta Visualization - PHATE 2/3...'); end
    S{end+1} = phate(sqrt(Xn), 't', 20, 'ndim', ndim, 'k', 15, 'pot_method', 'sqrt');
    
    if showwaitbar, gui.gui_waitbar_adv(fw, 5/nstep, 'Meta Visualization - PHATE 3/3...'); end
    S{end+1} = phate(sqrt(Xn), 't', 20, 'ndim', ndim, 'k', 30);
end

if showwaitbar, gui.gui_waitbar_adv(fw, 6/nstep, 'Meta Visualization - METAVIZ'); end
if usingmmfile
    [Y] = metaviz_memmap(S, ndim);
else
    [Y] = metaviz_tensor(S, ndim);
end
if showwaitbar, gui.gui_waitbar_adv(fw); end
end

%figure; scatter(Y(:,1),Y(:,2));

%{
disp('sTSNE2')
sTSNE2=tsne(Xn,Perplexity=15,NumDimensions=3);
disp('sTSNE3')
sTSNE3=tsne(Xn,Perplexity=50,NumDimensions=3);
disp('sTSNE4')
sTSNE4=tsne(data,Perplexity=30,NumDimensions=3);
disp('sTSNE5')
sTSNE5=tsne(data,Perplexity=15,NumDimensions=3);
disp('sTSNE6')
sTSNE6=tsne(data,Perplexity=50,NumDimensions=3);


disp('sUMAP2')
sUMAP2=run_umap_main(Xn,'n_components',3,'n_neighbors',30);
disp('sUMAP3')
sUMAP3=run_umap_main(Xn,'n_components',3,'n_neighbors',50);
disp('sUMAP4')
sUMAP4=run_umap_main(data,'n_components',3,'n_neighbors',15);
disp('sUMAP5')
sUMAP5=run_umap_main(data,'n_components',3,'n_neighbors',30);
disp('sUMAP6')
sUMAP6=run_umap_main(data,'n_components',3,'n_neighbors',50);

disp('sPHATE2')
sPHATE2 = phate(sqrt(Xn), 't', 20, 'ndim', 3, 'k', 30);
disp('sPHATE3')
sPHATE3 = phate(sqrt(Xn), 't', 20, 'ndim', 3, 'k', 50);
disp('sPHATE4')
sPHATE4 = phate(data, 't', 20, 'ndim', 3, 'k', 5);
disp('sPHATE5')
sPHATE5 = phate(data, 't', 20, 'ndim', 3, 'k', 30);
disp('sPHATE6')
sPHATE6 = phate(data, 't', 20, 'ndim', 3, 'k', 50);

%}

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

N = n * n * K;

mmf = tempname;
fileID = fopen(mmf, 'w');
for k = 1:K
    fwrite(fileID, zeros([n * n, 1]), 'single');
end
fclose(fileID);
m = memmapfile(mmf, 'Format', 'single', 'Writable', true);

% D=zeros(n,n,K,'single');
for k = 1:K
    d = pdist2(Sinput{k}, Sinput{k});
    %    D(:,:,k)=d./vecnorm(d);
    m.Data((n^2)*(k - 1)+1:(n^2)*(k)) = d ./ vecnorm(d);
end

%isequal(D(:),m.Data)

%%
m.Offset = 0;
for x = 1:n % n of cells
    d = reshape(m.Data(x:n:N), [n, K]);
    S = 1 - squareform(pdist(d', 'cosine'));

    %S1=1-squareform(pdist(squeeze(D(x,:,:))','cosine'));
    %isequal(S,S1)

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
        % isequal(s:t,D(:,l,k)')
        d = d + w(k, l)' .* m.Data(s:t);
        %d=d+w(k,l)'.*D(:,l,k);
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
for x = 1:n % n of cells
    S = 1 - squareform(pdist(squeeze(D(x, :, :))', 'cosine'));
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
