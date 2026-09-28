function [c, info] = e_snndpc(s, cluK, knnK)
%E_SNNDPC  Shared-nearest-neighbour density-peak clustering of an embedding.
%
%   C = PKG.E_SNNDPC(S, CLUK) partitions the rows of the embedding S into
%   CLUK clusters and returns a cluster index per cell.
%
%   C = PKG.E_SNNDPC(S, CLUK, KNNK) sets the neighbourhood size (default 4).
%
%   [C, INFO] = PKG.E_SNNDPC(...) also returns the intermediate quantities
%   RHO, DELTA, GAMMA, the shared-neighbour similarity DIST2 (sparse) and
%   the chosen CENTER indices, for comparison against the reference
%   implementation in EXTERNAL/ML_SNNDPC/SNNDPC_ORI.
%
%   This reproduces SnnDpc [DOI:10.1016/j.ins.2018.03.031] exactly and is
%   the function RUN.ML_SNNDPC now calls. SNNDPC_ORI is kept as the
%   reference the tests check this against; do not delete it.
%
%   WHY THIS EXISTS. The reference implementation is O(N^2) in time and
%   memory for a quantity that has O(N*K) nonzeros. On 8000 cells it spent
%   28.7 s of its 36.7 s in one scalar double loop over all N^2/2 pairs,
%   and held a 3.8 GB N-by-N CELL array of shared-neighbour index lists;
%   at 30000 cells that cell array alone would be ~54 GB. But DIST2(p,o)
%   is nonzero only when p and o are mutual K-nearest neighbours, and the
%   index lists were only ever used for their length and for one sum, so
%   the whole loop collapses into sparse products. With A the N-by-N kNN
%   indicator (K per row, the point itself included) and D the same
%   pattern holding the distances:
%
%       sharedCount = A*A'                    (exact, all pairs)
%       sum_q dist(p,q) + dist(o,q)   over q in KNN(p) n KNN(o)
%                   = (D*A' + A*D')(p,o)
%
%   which is where the 28.7 s went.
%
%   See also SC_SNNDPC, RUN.ML_SNNDPC, SC_CLUSTER_S.

if nargin < 3 || isempty(knnK), knnK = 4; end
if nargin < 2 || isempty(cluK), cluK = 10; end

N = size(s, 1);
K = knnK;
if K < 2
    error('pkg:e_snndpc:tooFewNeighbours', ...
        ['KNNK is %d. Use at least 2: a neighbourhood of size 1 holds ', ...
        'only the point itself.'], K);
end
if K + 1 > N
    error('pkg:e_snndpc:tooFewCells', ...
        ['KNNK is %d but there are only %d cells. KNNK+1 neighbours ', ...
        'must exist, so reduce KNNK below %d.'], K, N, N);
end

% Min-max normalise each coordinate, as the reference does. A constant
% column gives 0/0; the reference maps that to zero.
data = (s-min(s)) ./ (max(s) - min(s));
data(isnan(data)) = 0;

% ---- neighbourhoods -----------------------------------------------------
% Column 1 of KNNSEARCH is the point itself, matching the reference, where
% the K-neighbour set is the first K columns of the sorted distance order
% and so includes the point.
[nnIdx, nnD] = knnsearch(data, data, K=K+1);
rowIdx = repmat((1:N)', 1, K);
nbrIdx = nnIdx(:, 1:K);
A = sparse(rowIdx, nbrIdx, 1, N, N);
D = sparse(rowIdx, nbrIdx, nnD(:, 1:K), N, N);

sharedCount = A * A';
kDistNext = nnD(:, K+1);

% ---- shared-neighbour similarity ---------------------------------------
% Candidate pairs are the mutual neighbours: dist(p,o) < dist1Sort(p,K+1)
% already implies o is among the K nearest of p, so the reference's
% condition can only hold inside the pattern of A & A'.
T = D * A';
den = T + T';
[pIdx, oIdx] = find(triu(A&A', 1));
if isempty(pIdx)
    dist2 = sparse(N, N);
else
    lin = pIdx + (oIdx - 1) * N;
    cnt = full(sharedCount(lin));
    dsum = full(den(lin));
    dpo = vecnorm(data(pIdx, :)-data(oIdx, :), 2, 2);
    isMutual = dpo < min(kDistNext(pIdx), kDistNext(oIdx));
    v = zeros(numel(pIdx), 1);
    v(isMutual) = cnt(isMutual).^2 ./ dsum(isMutual);
    dist2 = sparse(pIdx, oIdx, v, N, N);
    dist2 = dist2 + dist2';
end

% ---- rho: sum of the K largest similarities in each row -----------------
rho = zeros(1, N);
[ri, ~, rv] = find(dist2);
if ~isempty(ri)
    [~, ord] = sortrows([ri, -rv]);
    ri = ri(ord);
    rv = rv(ord);
    perRow = accumarray(ri, 1, [N, 1]);
    rowStart = cumsum([0; perRow(1:end-1)]);
    rankInRow = (1:numel(ri))' - rowStart(ri);
    keep = rankInRow <= K;
    rho = accumarray(ri(keep), rv(keep), [N, 1])';
end

% ---- delta: distance to the nearest point of higher rho -----------------
% Inherently all-pairs, but done in chunks against the coordinates rather
% than against a materialised N-by-N distance matrix.
kDistSum = sum(nnD(:, 1:K), 2);
[~, rhoOrder] = sort(rho, 'descend');
ordData = data(rhoOrder, :);
ordSum = kDistSum(rhoOrder);
deltaOrd = inf(N, 1);
selectOrd = zeros(N, 1);
chunk = max(1, floor(4e6/N));
for a = 1:chunk:N
    b = min(a+chunk-1, N);
    if b < 2, continue; end
    dChunk = pdist2(ordData(a:b, :), ordData(1:b, :));
    cand = dChunk .* (ordSum(a:b) + ordSum(1:b)');
    % Only points of strictly higher rho are candidates.
    cand((1:b) >= (a:b)') = inf;
    [mv, mi] = min(cand, [], 2);
    deltaOrd(a:b) = mv;
    selectOrd(a:b) = mi;
end
selectOrd(1) = 0;   % the highest-rho point has no candidate
deltaOrd(1) = -1;
deltaOrd(1) = max(deltaOrd);
delta = zeros(1, N);
delta(rhoOrder) = deltaOrd;
deltaSelect = zeros(1, N);
nz = selectOrd > 0;
deltaSelect(rhoOrder(nz)) = rhoOrder(selectOrd(nz));

gamma = rho .* delta;

% ---- pick centres -------------------------------------------------------
gammaSort = sort(gamma, 'descend');
yCut = mean(gammaSort(cluK:cluK+1));
center = find(gamma > yCut);
NC = numel(center);
cluster = -ones(1, N);
cluster(center) = 1:NC;

% ---- assign the inevitable subordinate points ---------------------------
% Same FIFO order as the reference, which offered centres at the head of a
% java LinkedList and polled from the tail. A preallocated array does it
% without the java round trip.
scNbr = full(sharedCount(sub2ind([N, N], ...
    repmat((1:N)', 1, K-1), nnIdx(:, 2:K))));
queue = zeros(1, N);
queue(1:NC) = center;
head = 1;
tail = NC;
while head <= tail
    this = queue(head);
    head = head + 1;
    nbr = nnIdx(this, 2:K);
    take = cluster(nbr) < 0 & scNbr(this, :) >= K/2;
    for next = nbr(take)
        % Re-test: an earlier NEXT in this same row may have claimed it.
        if cluster(next) < 0
            cluster(next) = cluster(this);
            tail = tail + 1;
            queue(tail) = next;
        end
    end
end

% ---- assign the possible subordinate points -----------------------------
kCur = K;
unas = find(cluster < 0);
while ~isempty(unas)
    nb = nnIdx(unas, 2:kCur);
    cl = cluster(nb);
    valid = cl > 0;
    if any(valid, 'all')
        rIdx = repmat((1:numel(unas))', 1, kCur-1);
        rSub = rIdx(valid);
        cSub = cl(valid);
        recog = accumarray([rSub(:), cSub(:)], 1, [numel(unas), NC]);
        whichK = max(recog(:));
    else
        whichK = 0;
    end
    if whichK > 0
        % Column-major FIND order, so a point tied between clusters takes
        % the highest cluster index -- what the reference did.
        [whichPoint, whichCluster] = find(recog == whichK);
        cluster(unas(whichPoint)) = whichCluster;
        unas = find(cluster < 0);
    else
        kCur = kCur + 1;
        if kCur > N
            error('pkg:e_snndpc:unassignable', ...
                ['%d cells could not be attached to any cluster. The ', ...
                'embedding may hold duplicate or disconnected points.'], ...
                numel(unas));
        end
        if kCur > size(nnIdx, 2)
            nq = min(N, 2*size(nnIdx, 2));
            nnIdx = knnsearch(data, data, K=nq);
        end
    end
end

c = cluster(:);

if nargout > 1
    info = struct('NC', NC, 'K', K, 'dist2', dist2, 'rho', rho, ...
        'delta', delta, 'deltaSelect', deltaSelect, 'gamma', gamma, ...
        'cluster', c, 'center', center(:), 'sharedCount', sharedCount);
end
end
