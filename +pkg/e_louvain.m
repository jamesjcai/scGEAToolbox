function [c, Q] = e_louvain(A, gamma, maxPasses)
%E_LOUVAIN  Louvain community detection on a weighted graph.
%
%   C = PKG.E_LOUVAIN(A) partitions the nodes of the sparse symmetric
%   adjacency matrix A into communities that greedily maximise modularity,
%   and returns a community index per node.
%
%   C = PKG.E_LOUVAIN(A, GAMMA) sets the resolution (default 1). Larger
%   values give more, smaller communities; the number of communities is
%   not set directly. SC_LOUVAIN tunes GAMMA when a target count is wanted.
%
%   [C, Q] = PKG.E_LOUVAIN(...) also returns the modularity of the result,
%   at the same resolution.
%
%   A may carry self-loops, which the aggregated levels of the algorithm
%   produce; they count towards a node's degree but never towards the
%   links it has into another community.
%
%   Nodes are visited in index order rather than a random one, so repeated
%   calls on the same graph return the same partition.
%
%   Reference: Blondel et al. (2008), J Stat Mech P10008,
%   DOI:10.1088/1742-5468/2008/10/P10008. The resolution enters as in
%   Reichardt & Bornholdt (2006), DOI:10.1103/PhysRevE.74.016110.
%
%   See also SC_LOUVAIN, SC_KNNGRAPH, SC_CLUSTER_S.

if nargin < 3 || isempty(maxPasses), maxPasses = 100; end
if nargin < 2 || isempty(gamma), gamma = 1; end

A = double(A);
if ~issparse(A), A = sparse(A); end
n0 = size(A, 1);
if size(A, 2) ~= n0
    error('pkg:e_louvain:notSquare', ...
        'A is %d-by-%d. The adjacency matrix must be square.', ...
        size(A, 1), size(A, 2));
end
A = (A + A') / 2;

twoM = full(sum(A, 'all'));
if twoM <= 0
    % No edges, so no partition is better than any other.
    c = (1:n0)';
    Q = 0;
    return
end

map = (1:n0)';
W = A;

while true
    n = size(W, 1);
    k = full(sum(W, 2));
    Wnd = W - diag(sparse(diag(W)));

    % Flatten the sparse structure once. The inner loop below runs once
    % per node per pass, and repeated sparse column indexing there costs
    % more than the moves themselves.
    [nbrAll, colAll, wAll] = find(Wnd);
    % FIND returns rows, not columns, for a 1-by-1 matrix, which is what W
    % becomes once every node has merged into one community.
    nbrAll = nbrAll(:);
    colAll = colAll(:);
    wAll = wAll(:);
    perCol = accumarray(colAll, 1, [n, 1]);
    colStart = cumsum([1; perCol]);

    comm = (1:n)';
    sigTot = k;
    wTo = zeros(n, 1);
    seen = false(n, 1);
    touched = zeros(n, 1);

    for pass = 1:maxPasses
        moved = 0;
        for i = 1:n
            lo = colStart(i);
            hi = colStart(i+1) - 1;
            if hi < lo, continue; end

            ci = comm(i);
            ki = k(i);
            sigTot(ci) = sigTot(ci) - ki;

            % Weight from i into each neighbouring community.
            nTouched = 0;
            for t = lo:hi
                ct = comm(nbrAll(t));
                if ~seen(ct)
                    seen(ct) = true;
                    nTouched = nTouched + 1;
                    touched(nTouched) = ct;
                end
                wTo(ct) = wTo(ct) + wAll(t);
            end

            % Staying put is always a candidate, and wins ties.
            bestGain = wTo(ci) - gamma*ki*sigTot(ci)/twoM;
            cBest = ci;
            for t = 1:nTouched
                ct = touched(t);
                g = wTo(ct) - gamma*ki*sigTot(ct)/twoM;
                if g > bestGain
                    bestGain = g;
                    cBest = ct;
                end
            end

            for t = 1:nTouched
                ct = touched(t);
                wTo(ct) = 0;
                seen(ct) = false;
            end

            comm(i) = cBest;
            sigTot(cBest) = sigTot(cBest) + ki;
            if cBest ~= ci, moved = moved + 1; end
        end
        if moved == 0, break; end
    end

    [~, ~, newComm] = unique(comm);
    nc = max(newComm);
    map = newComm(map);
    if nc == n, break; end

    P = sparse(1:n, newComm, 1, n, nc);
    W = P' * W * P;
end

c = map;

if nargout > 1
    P = sparse(1:n0, c, 1, n0, max(c));
    B = P' * A * P;
    kc = full(sum(B, 2));
    Q = sum(full(diag(B))/twoM-gamma*(kc / twoM).^2);
end
end
