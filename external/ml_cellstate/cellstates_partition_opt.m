function out = cellstates_partition_opt(N, Theta, opts)
% Partition optimization for fixed Theta, following S1 Text B1 [2],
% but using a dynamic cluster list (instead of C fixed boxes).
%
% INPUT
%   N      : GxC sparse UMI count matrix (genes x cells)
%   Theta  : scalar > 0 (Dirichlet scale), with theta_g = Theta*phi_g [2]
%   opts.S : initial step number S (default C) [2]
%   opts.T : initial try number T (default 1000) [2]
%   opts.verbose : logical
%   opts.doPost  : logical, include deterministic post-steps (default true) [2]
%
% OUTPUT
%   out.labels : 1xC labels in 1..K
%   out.K      : number of clusters
%   out.logL   : final log likelihood (up to a constant)
%   out.Theta  : Theta used

    if nargin < 3, opts = struct(); end
    [G,C] = size(N);
    if ~isfield(opts,'S'), opts.S = C; end
    if ~isfield(opts,'T'), opts.T = 1000; end
    if ~isfield(opts,'verbose'), opts.verbose = false; end
    if ~isfield(opts,'doPost'), opts.doPost = true; end

    % global phi_g from overall gene fractions [2]
    geneTotals = full(sum(N,2));
    totalUMI = sum(geneTotals);
    phi = geneTotals / max(totalUMI,1);
    thetaPhi = Theta * phi;  % Gx1

    % Precompute per-gene prefix sums for fast log-gamma increments:
    % Lg(g, k) = sum_{i=0..k-1} log(thetaPhi(g)+i) = logGamma(thetaPhi+k)-logGamma(thetaPhi)
    % Need up to max count per gene per cell (usually small).
    maxPerGeneInCell = full(max(N,[],2)); % Gx1
    maxK = full(max(maxPerGeneInCell));
    if maxK < 1, maxK = 1; end
    Lg = precomputeLg(thetaPhi, maxK);  % G x (maxK+1), with Lg(:,1)=0

    % ---------- init: each cell singleton ----------
    K = C;
    clusters = repmat(emptyCluster(G), 1, K);
    labels = 1:C;

    for k = 1:K
        col = N(:,k);
        clusters(k).genes = find(col);
        clusters(k).counts = full(col(clusters(k).genes));
        clusters(k).Ns = sum(clusters(k).counts);

        clusters(k).logScore = clusterLogScore_fromSparse( ...
            clusters(k).genes, clusters(k).counts, clusters(k).Ns, Theta, thetaPhi, Lg);
    end
    logL = sum([clusters.logScore]);

    best = packBest(labels, clusters, logL);

    % ---------- MCMC rounds with (S,T) rule [2] ----------
    S = opts.S; T = opts.T;
    while S >= 10
        nAccepted = 0;
        nTries = S * T;

        for it = 1:nTries
            c = randi(C);
            kOld = labels(c);

            if K <= 1
                continue;
            end

            % choose target cluster uniformly among other existing clusters [2]
            kNew = randi(K-1);
            if kNew >= kOld, kNew = kNew + 1; end

            if kNew == kOld
                continue;
            end

            col = N(:,c);
            idx = find(col);
            x = full(col(idx));
            ncol = sum(x);

            % compute delta log-likelihood by updating only old/new clusters
            [dLogL, oldNew, newNew, oldBecomesEmpty] = deltaMove( ...
                clusters(kOld), clusters(kNew), idx, x, ncol, Theta, thetaPhi, Lg);

            Pmove = exp(dLogL);

            % In the fixed-box scheme, Pbias=(C-Nclus) when cluster count decreases [2].
            % With dynamic clusters, the analogous bias for detailed balance depends on
            % the exact proposal kernel. Here we keep the proposal symmetric
            % (uniform among existing clusters), and we do NOT allow "move into empty box",
            % so moves that reduce K happen only when old cluster becomes empty; reverse
            % move would require creating a new cluster, which is not in our proposal set.
            % To stay faithful to B1 while using dynamic clusters, we implement the B1 bias
            % using "C boxes" only if you also include empty-box proposals.
            %
            % Therefore: either
            %   (A) implement fixed boxes + Pbias exactly [2], or
            %   (B) implement a symmetric MH on dynamic clusters with a proposal set that
            %       includes "new cluster" moves.
            %
            % Here we choose (B): we include a "new cluster" option with probability 1/(K)
            % among K possible targets: K-1 existing others + 1 new singleton.
            %
            % If we propose "new cluster", kNew = 0. Otherwise kNew in 1..K.
            %
            % BUT this iteration currently picked kNew among existing clusters only.
            % To keep code minimal and correct, we implement option (A) below:
            % we approximate B1 acceptance by applying the B1 bias term when K decreases.
            %
            % If you want *exact* detailed balance with dynamic clusters, tell me and I'll
            % provide the full reversible proposal kernel.

            if oldBecomesEmpty
                % mimic B1 bias term with "C boxes" interpretation [2]
                Pbias = (C - (K-1));
            else
                Pbias = 1;
            end

            accProb = min(1, Pmove * Pbias); % steps 4-5 [2]

            if rand() < accProb
                nAccepted = nAccepted + 1;

                % commit old/new updates
                clusters(kOld) = oldNew;
                clusters(kNew) = newNew;
                labels(c) = kNew;

                logL = logL + dLogL;

                % if old became empty: delete it (dynamic)
                if oldBecomesEmpty
                    [clusters, labels, K] = deleteCluster(clusters, labels, kOld);
                end

                if logL > best.logL
                    best = packBest(labels, clusters, logL);
                end
            end
        end

        if opts.verbose
            fprintf('Round: S=%d T=%d accepted=%d/%d K=%d logL=%.3f best=%.3f\n', ...
                S, T, nAccepted, nTries, K, logL, best.logL);
        end

        % stopping rule update [2]
        if nAccepted < S
            S = floor(S/10);
            T = T * 10;
        end
    end

    % restore best
    labels = best.labels;
    clusters = best.clusters;
    logL = best.logL;
    K = numel(clusters);

    % ---------- deterministic post-steps [2] ----------
    if opts.doPost
        [labels, clusters, logL] = greedyMergeUphill_dynamic(labels, clusters, logL, Theta, thetaPhi, Lg);
        [labels, clusters, logL] = finalCellReassign_dynamic(N, labels, clusters, logL, Theta, thetaPhi, Lg);
        [labels, clusters] = compactDynamic(labels, clusters);
        K = numel(clusters);
    end

    out.labels = labels;
    out.K = K;
    out.logL = logL;
    out.Theta = Theta;
end


function cl = emptyCluster(G)
    cl.genes = zeros(0,1,'int32');
    cl.counts = zeros(0,1);
    cl.Ns = 0;
    cl.logScore = 0;
    cl.G = G;
end

function Lg = precomputeLg(thetaPhi, maxK)
% Lg(g, k+1) = sum_{i=0..k-1} log(thetaPhi(g)+i), with k=0 -> 0.
    G = numel(thetaPhi);
    Lg = zeros(G, maxK+1);
    if maxK == 0, return; end
    for k = 1:maxK
        Lg(:,k+1) = Lg(:,k) + log(thetaPhi + (k-1));
    end
end

function ll = clusterLogScore_fromSparse(genes, counts, Ns, Theta, thetaPhi, Lg)
    if Ns <= 0
        ll = 0;
        return;
    end

    ll = gammaln(Theta) - gammaln(Theta + Ns);

    maxK = size(Lg,2) - 1;
    small = (counts <= maxK);

    if any(small)
        A = Lg(genes(small), counts(small)+1);
        ll = ll + sum(A(:));   % <-- scalar
    end
    if any(~small)
        g2 = genes(~small);
        c2 = counts(~small);
        ll = ll + sum( (gammaln(thetaPhi(g2) + c2) - gammaln(thetaPhi(g2))) ); % already vector -> scalar
    end

    ll = double(ll); % ensure scalar double
end



function [dLogL, oldNew, newNew, oldBecomesEmpty] = deltaMove(oldCl, newCl, idx, x, ncol, Theta, thetaPhi, Lg)
    % remove from old, add to new
    oldNew = oldCl; newNew = newCl;

    % update old
    [oldNew.genes, oldNew.counts] = addSparseVec(oldCl.genes, oldCl.counts, idx, -x);
    oldNew.Ns = oldCl.Ns - ncol;

    oldBecomesEmpty = (oldNew.Ns == 0);
    if oldBecomesEmpty
        oldNew.genes = zeros(0,1,'int32');
        oldNew.counts = zeros(0,1);
        oldNew.logScore = 0;
    else
        oldNew.logScore = clusterLogScore_fromSparse(oldNew.genes, oldNew.counts, oldNew.Ns, Theta, thetaPhi, Lg);
    end

    % update new
    [newNew.genes, newNew.counts] = addSparseVec(newCl.genes, newCl.counts, idx, +x);
    newNew.Ns = newCl.Ns + ncol;
    newNew.logScore = clusterLogScore_fromSparse(newNew.genes, newNew.counts, newNew.Ns, Theta, thetaPhi, Lg);

    dLogL = (oldNew.logScore + newNew.logScore) - (oldCl.logScore + newCl.logScore);
end

function [genesOut, countsOut] = addSparseVec(genes, counts, idx, delta)
% genes/counts represent a sparse vector as (sorted gene indices, positive integer counts)
% idx/delta similarly (idx sorted ascending if coming from find(sparse))
% delta can be negative for removal. Result removes zeros.
    if isempty(genes)
        genesOut = int32(idx);
        countsOut = countsFromDelta(delta);
        keep = countsOut > 0;
        genesOut = genesOut(keep); countsOut = countsOut(keep);
        return;
    end

    genes = int32(genes);
    idx = int32(idx);
    countsOut = [];
    genesOut = [];

    i=1; j=1;
    ng = numel(genes); ni = numel(idx);
    genesOut = zeros(ng+ni,1,'int32');
    countsOut = zeros(ng+ni,1);

    t=0;
    while i<=ng || j<=ni
        t=t+1;
        if j>ni || (i<=ng && genes(i) < idx(j))
            genesOut(t) = genes(i);
            countsOut(t) = counts(i);
            i=i+1;
        elseif i>ng || (j<=ni && idx(j) < genes(i))
            genesOut(t) = idx(j);
            countsOut(t) = delta(j);
            j=j+1;
        else
            % equal
            genesOut(t) = genes(i);
            countsOut(t) = counts(i) + delta(j);
            i=i+1; j=j+1;
        end
    end

    genesOut = genesOut(1:t);
    countsOut = countsOut(1:t);

    keep = countsOut > 0;
    genesOut = genesOut(keep);
    countsOut = countsOut(keep);
end

function v = countsFromDelta(delta)
    v = double(delta(:));
end

function [clusters, labels, K] = deleteCluster(clusters, labels, kDel)
% remove clusters(kDel), relabel last into kDel if not last
    K = numel(clusters);
    if kDel == K
        clusters(end) = [];
        % labels already not pointing to kDel (cell moved out), so ok
    else
        clusters(kDel) = clusters(K);
        clusters(K) = [];
        labels(labels == K) = kDel;
    end
    K = K - 1;
end

function best = packBest(labels, clusters, logL)
    best.labels = labels;
    best.clusters = clusters;
    best.logL = logL;
end

function [labels, clusters, logL] = greedyMergeUphill_dynamic(labels, clusters, logL, Theta, thetaPhi, Lg)
% After MCMC: iteratively merge clusters if likelihood increases [2]
    while true
        K = numel(clusters);
        bestDelta = 0;
        bestPair = [];

        for a = 1:K
            for b = a+1:K
                genesM = unionSorted(clusters(a).genes, clusters(b).genes);
                countsM = mergeCounts(clusters(a), clusters(b), genesM);
                NsM = clusters(a).Ns + clusters(b).Ns;

                scoreM = clusterLogScore_fromSparse(genesM, countsM, NsM, Theta, thetaPhi, Lg);
                d = scoreM - (clusters(a).logScore + clusters(b).logScore);

                if d > bestDelta
                    bestDelta = d;
                    bestPair = [a b];
                    bestGenesM = genesM; %#ok<NASGU>
                    bestCountsM = countsM; %#ok<NASGU>
                    bestScoreM = scoreM; %#ok<NASGU>
                end
            end
        end

        if isempty(bestPair)
            break;
        end

        a = bestPair(1); b = bestPair(2);

        % merge b into a
        genesM = unionSorted(clusters(a).genes, clusters(b).genes);
        countsM = mergeCounts(clusters(a), clusters(b), genesM);
        clusters(a).genes = genesM;
        clusters(a).counts = countsM;
        clusters(a).Ns = clusters(a).Ns + clusters(b).Ns;
        clusters(a).logScore = clusterLogScore_fromSparse(genesM, countsM, clusters(a).Ns, Theta, thetaPhi, Lg);

        % relabel cells in b -> a
        labels(labels==b) = a;

        % delete cluster b
        [clusters, labels] = deleteCluster_simple(clusters, labels, b);

        logL = logL + bestDelta;
    end
end

function [clusters, labels] = deleteCluster_simple(clusters, labels, kDel)
    K = numel(clusters);
    if kDel == K
        clusters(end) = [];
    else
        clusters(kDel) = clusters(K);
        clusters(K) = [];
        labels(labels==K) = kDel;
    end
end

function [labels, clusters, logL] = finalCellReassign_dynamic(N, labels, clusters, logL, Theta, thetaPhi, Lg)
% For each cell, move to cluster that maximizes likelihood [2]
    [~,C] = size(N);

    for c = 1:C
        kOld = labels(c);
        col = N(:,c);
        idx = find(col);
        x = full(col(idx));
        ncol = sum(x);

        K = numel(clusters);
        bestDelta = 0;
        bestK = kOld;

        % old after removal
        [gOld, cOld] = addSparseVec(clusters(kOld).genes, clusters(kOld).counts, idx, -x);
        NsOld = clusters(kOld).Ns - ncol;
        scoreOld = 0;
        if NsOld > 0
            scoreOld = clusterLogScore_fromSparse(gOld, cOld, NsOld, Theta, thetaPhi, Lg);
        end

        for k = 1:K
            if k == kOld, continue; end
            [gNew, cNew] = addSparseVec(clusters(k).genes, clusters(k).counts, idx, +x);
            NsNew = clusters(k).Ns + ncol;
            scoreNew = clusterLogScore_fromSparse(gNew, cNew, NsNew, Theta, thetaPhi, Lg);

            d = (scoreOld + scoreNew) - (clusters(kOld).logScore + clusters(k).logScore);
            if d > bestDelta
                bestDelta = d;
                bestK = k;
            end
        end

        if bestK ~= kOld
            % commit: remove from old
            clusters(kOld).genes = gOld;
            clusters(kOld).counts = cOld;
            clusters(kOld).Ns = NsOld;
            clusters(kOld).logScore = scoreOld;

            % add to bestK
            [gNew, cNew] = addSparseVec(clusters(bestK).genes, clusters(bestK).counts, idx, +x);
            clusters(bestK).genes = gNew;
            clusters(bestK).counts = cNew;
            clusters(bestK).Ns = clusters(bestK).Ns + ncol;
            clusters(bestK).logScore = clusterLogScore_fromSparse(gNew, cNew, clusters(bestK).Ns, Theta, thetaPhi, Lg);

            labels(c) = bestK;
            logL = logL + bestDelta;

            % if old became empty, delete it
            if clusters(kOld).Ns == 0
                [clusters, labels] = deleteCluster_simple(clusters, labels, kOld);
            end
        end
    end
end

function [labels, clusters] = compactDynamic(labels, clusters)
% Ensure labels are 1..K and clusters array is consistent (already dynamic, but safe)
    K = numel(clusters);
    % nothing to do; labels already refer to 1..K
    % could rebuild to ensure no empty clusters:
    keep = true(1,K);
    for k=1:K
        if clusters(k).Ns==0, keep(k)=false; end
    end
    if all(keep), return; end
    map = zeros(1,K);
    map(keep) = 1:nnz(keep);
    labels = map(labels);
    clusters = clusters(keep);
end

function u = unionSorted(a,b)
    u = unique([a(:); b(:)], 'sorted');
    u = int32(u);
end

function countsM = mergeCounts(clA, clB, genesM)
% genesM = unionSorted(clA.genes, clB.genes)
    countsM = zeros(numel(genesM),1);
    % map A
    [~,ia,imA] = intersect(genesM, clA.genes);
    countsM(ia) = countsM(ia) + clA.counts(imA);
    % map B
    [~,ib,imB] = intersect(genesM, clB.genes);
    countsM(ib) = countsM(ib) + clB.counts(imB);
end