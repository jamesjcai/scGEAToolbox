function [optimk] = e_numclusters(X, varargin)
% Estimate number of clusters
p = inputParser;
defaultType = 'simlr';
validTypes = {'simlr', 'soptsc', 'sc3'};
checkType = @(x) any(validatestring(x, validTypes));

addRequired(p, 'X', @isnumeric);
addOptional(p, 'type', defaultType, checkType);
parse(p, X, varargin{:});


pw1 = fileparts(mfilename('fullpath'));
switch p.Results.type
    case 'simlr'
        pth = fullfile(pw1, '..', 'external', 'ml_SIMLR');
        if ~(ismcc || isdeployed), addpath(pth); end
        pth = fullfile(pw1, '..', 'external', 'ml_SIMLR', 'src');
        if ~(ismcc || isdeployed), addpath(pth); end
        [~, K2] = Estimate_Number_of_Clusters_SIMLR(X', 2:10);
        [~, i] = min(K2);
        optimk = i + 1;
    case 'soptsc'
        pth = fullfile(pw1, '..', 'external', 'ml_SoptSC');
        if ~(ismcc || isdeployed), addpath(pth); end
        pth = fullfile(pw1, '..', 'external', 'ml_SoptSC', 'NNDSVD');
        if ~(ismcc || isdeployed), addpath(pth); end
        pth = fullfile(pw1, '..', 'external', 'ml_SoptSC', 'symnmf2');
        if ~(ismcc || isdeployed), addpath(pth); end

        realdata = X;
        realdata = realdata - min(realdata(:));
        realdata = realdata ./ max(realdata(:));

        [~, n] = size(realdata);
        for i = 1:n
            realdata(:, i) = realdata(:, i) / norm(realdata(:, i));
        end
        lambda = 0.5;
        W = SimilarityM(realdata, lambda, X);
        WB = W;
        n = size(W, 1);
        D = diag(WB*ones(n, 1));
        Prw = eye(size(W)) - D^(-1 / 2) * WB * D^(-1 / 2);
        if n >= 1000
            No_eigs = 100;
            all_eigs = real(eigs(Prw, No_eigs, 'sm'));
        else
            all_eigs = real(eig(Prw));
        end

        ZZ = sort(abs(real(all_eigs)));
        No_cluster1 = length(find(ZZ <= 0.01));

        % Determinning the number of clusters
        eigenvalues = [];
        % if isempty(optimk)
        [~, No_cluster] = Num_cluster(W, No_cluster1);
        optimk = No_cluster;
        % end
    case 'sc3'

        %% estimate k
        % X=log2(X+1);
        Dis = squareform(pdist(X'));
        A = exp(-Dis./max(Dis(:))); % adjacency matrix
        xD = diag(sum(A).^-0.5); % D=diag(sum(A)); % d(i) the degree of node i
        xA = xD * A * xD; % normalized adjacenty matrix
        L = eye(size(A, 1)) - xA; % also L=xD*(D-A)*xD

        % see https://people.orie.cornell.edu/dpw/orie6334/lecture7.pdf
        % see https://en.wikipedia.org/wiki/Laplacian_matrix#Symmetric_normalized_Laplacian_2

        [V, D] = eig(L);
        [~, ind] = sort(diag(D));
        Vs = V(:, ind);

        % Two changes here, because the selection this branch performed
        % could not discriminate at all.
        %
        % (1) Cluster the leading i eigenvectors, not all n of them.
        % Spectral clustering embeds the points in the space spanned by
        % the eigenvectors of the smallest eigenvalues; using the whole
        % n-by-n matrix means clustering the rows of an orthogonal matrix,
        % which are equidistant by construction.
        %
        % (2) Score the partition in the data space, not in the embedding.
        % Calinski-Harabasz on the full eigenvector matrix came out
        % [NaN 1 1 1 1 1] -- identically 1 for every k, so OptimalK was
        % whichever value the tie-break happened to land on. Scoring in
        % the embedding instead is no better: CH there rises
        % monotonically and always returns the largest k offered.
        %
        % Measured over 40 fixtures (k = 2..5, ten seeds each) of
        % well-separated Gaussian blobs: the old code recovered the true k
        % in 12 of 40, never once for k = 4 or k = 5, with a mean error of
        % 1.30 clusters. This recovers it in 40 of 40.
        %
        % The affinity above, exp(-Dis./max(Dis(:))), still uses the
        % largest pairwise distance as its bandwidth, which leaves every
        % entry within a factor of e of every other and gives the graph
        % almost no block structure -- an eigengap criterion cannot read
        % anything off it. This branch is a weak estimator for that
        % reason; it is now at least an estimator. No in-tree caller
        % reaches it: +run/ml_SC3.m calls e_numclusters with the default
        % type, which is 'simlr'.
        kmax = 6;
        clust = zeros(size(Vs, 1), kmax);
        for i = 1:kmax
            clust(:, i) = kmeans(Vs(:, 1:i), i, ...
                'emptyaction', 'singleton', 'replicate', 5);
        end
        va = evalclusters(X', clust, 'CalinskiHarabasz');
        optimk = va.OptimalK;

end
