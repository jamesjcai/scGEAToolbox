function [aln0, aln1] = i_ma(A0, A1, ndim, laplacian, symmetrize, opts)
% MA - manifold alignment
%
% [aln0, aln1] = ten.i_ma(A0, A1)
% [aln0, aln1] = ten.i_ma(A0, A1, ndim, laplacian)
% [aln0, aln1] = ten.i_ma(A0, A1, ndim, laplacian, symmetrize)
% [aln0, aln1] = ten.i_ma(A0, [], ndim, laplacian, true, Knockouts=K)
%
% KNOCKOUTS=K aligns A0 with each of numel(K) knocked-out copies of itself,
% copy b being A0 with row K(b) zeroed, and returns ALN0 and ALN1 as
% numel(K)-by-1 cells. It is ten.i_knk's alignment for many genes at once:
% the same result as ten.i_ma(A0, A1_b, ...) for each, but A0's part of
% the work is done once and each knockout's change is applied as the
% sparse row-and-column update it is (ten.knknullcontrol at 3000 genes,
% 60 controls: 23 s -> 13 s). Symmetrized alignment only.
%
% LAPLACIAN is "normalized" (default) or "unnormalized".
%
% SYMMETRIZE (default true) symmetrizes the joint adjacency W; see below for
% why. FALSE reproduces manifoldAlignment() of the R packages, which works on
% W as it is: with "unnormalized" the result then matches R up to the sign of
% each column (which drdist does not see). Only asymmetric input is affected,
% i.e. scTenifoldKnk, whose knocked-out row makes A1 asymmetric. The reference
% mode of TEN.SCTENIFOLDKNK uses it; the Laplacian is then not symmetric, so
% eigs takes the non-symmetric path and the real parts of its modes are kept,
% as R does with Re().
%
% THE DEFAULT CHANGED TO THE NORMALIZED LAPLACIAN (2026-09-24). This used the
% unnormalized L = D - W, as the reference scTenifoldNet does. Its smallest
% modes concentrate on weakly connected nodes, so drdist was driven by
% peripheral genes whose few surviving edges are noise, not by the genes whose
% wiring changed. The symmetric normalized L = I - D^(-1/2) W D^(-1/2)
% divides that degree effect out.
%
% Measured with the full pipeline (10 bootstrap networks of 500 cells, CP
% denoising, 1000 genes, 3 seeds), on real Alpha cells split into two halves
% with 20 genes shuffled across cells in one half:
%   planted-gene AUROC, PCR networks    0.92 -> 0.93 (full shuffle),
%                                       0.66 -> 0.69 (30% of cells)
%   genes at FDR < 0.05 on a null split (no change) 12 -> 5
% On three real two-condition comparisons (ovarian cancer cells and
% fibroblasts, co-culture vs monoculture; tumour vs normal) the ranking
% correlates 0.93-0.96 (Spearman) with the unnormalized one, while null-split
% hits fall from 10-18 to 0-0.7. The chi2 p-values in ten.i_dr get stricter
% under the normalized Laplacian: expect fewer genes at FDR < 0.05.
%
% The same benchmark is why net.idsnet is not scTenifoldNet's network builder:
% it builds a network 20-30x faster than PCR, but even with this Laplacian it
% scored 0.79 against PCR's 0.93 on the planted genes, and on the real data
% only 14-32% of its top 50 genes matched PCR's.
%
% scTenifoldKnk (ten.i_knk, ndim = 2) uses the new default too. On its own
% denoised networks (1000 genes, three real populations, 40 random knockouts
% each) the KO gene ranked first in every knockout (85-100% before), the
% top 50 held 84-96% of the 50 genes losing the most edge weight (53-64%
% before), and two knockouts' top-50 lists overlapped 26-37% (51-59%
% before). The price is the tail: across all genes, drdist tracks the
% removed edge weight at Spearman 0.53-0.59, against 0.92-0.93 before, so a
% rank-based analysis over the whole Knk list is better served by
% "unnormalized". Re-measured on 2026-09-25 with the corrected (outgoing)
% knockout of TEN.I_KNK: the same figures to within a few percent.
%
% SYMMETRIZE THE ADJACENCY BEFORE USING "normalized". With SYMMETRIZE =
% false the normalized Laplacian loses the knockout: on R-pipeline Knk
% networks its top 50 held 6-11% of the strongest removed edges, against 86%
% on the same Alpha-cell network symmetrized.
%
% The modes returned are eigenvectors of the normalized L itself, not the
% D^(-1/2)-rescaled random-walk ones; the benchmark above measured these.
%
% Pass "unnormalized" to reproduce scTenifoldNet output from before this date.

arguments
    A0 {mustBeNumeric}
    A1 {mustBeNumeric}
    ndim (1, 1) double {mustBeInteger, mustBePositive} = 30
    laplacian (1, 1) string = "normalized"
    symmetrize (1, 1) logical = true
    opts.Knockouts double {mustBeInteger, mustBePositive} = []
end
if ~any(laplacian == ["normalized", "unnormalized"])
    error('ten:i_ma:badLaplacian', ...
        'LAPLACIAN must be "normalized" or "unnormalized". Received "%s".', ...
        laplacian);
end
if ~isempty(opts.Knockouts)
    if ~symmetrize || ~isempty(A1)
        error('ten:i_ma:knockoutsUsage', ...
            'KNOCKOUTS takes A1 = [] and the symmetrized alignment.');
    end
    [aln0, aln1] = i_knockouts(A0, opts.Knockouts, ndim, laplacian);
    return
end
mu = 0.9;

W1 = A0 + 1;
W2 = A1 + 1;

% W12 = eye(p1, p2) couples gene i of one network with gene i of the other.
mu = mu * (sum(W1(:)) + sum(W2(:))) / (2 * min(size(W1, 2), size(W2, 2)));

% SYMMETRIZE THE ADJACENCY, NOT THE LAPLACIAN.
%
% W is asymmetric whenever A1 is, and under scTenifoldKnk A1 is asymmetric BY
% CONSTRUCTION: ten.i_knk zeroes the knocked-out gene's row, which leaves its
% column intact. That happens even when A0 is perfectly symmetric - measured
% max|A1-A1'| = 0.153 starting from an exactly symmetric A0. scTenifoldNet
% never sees it, because it symmetrizes both matrices before calling here.
%
% This used to read
%
%     D = sum(abs(W));
%     L = diag(D) - W;
%     L = 0.5 * (L + L.');    % symmetric by construction; make it exactly so
%
% and the comment was wrong on both counts on the Knk path. The column sum on
% the first line is right and matches the reference (R's
% diag(W) <- -apply(W,2,sum)): with W >= 0 it makes the constant vector the
% exact LEFT null vector, so L has an exact zero eigenvalue even while
% asymmetric. Symmetrizing the Laplacian afterwards destroyed that. The result
% was L_proper + diag(delta), delta_j = (colsum_j - rowsum_j)/2, which has
% neither row nor column sums zero - the Laplacian of no graph. Because
% 1'L1 = 0 while L1 =/= 0, such a matrix is strictly indefinite: the measured
% smallest eigenvalue was -1.4e-05 where a real Laplacian gives ~1e-12.
%
% It changed answers. On a 5005-gene network the knocked-out gene came out at
% rank 3, behind two high-degree hub genes; symmetrizing W instead puts it at
% rank 1 where it belongs, and the two rankings correlate at Spearman -0.23.
%
% Symmetrizing W keeps what the old line was actually buying - a symmetric L,
% so eigs takes the Lanczos path and returns real vectors rather than the
% complex pairs an asymmetric L can produce - while also giving a genuine PSD
% Laplacian with L*1 = 0 and 1'L = 0. This is what ten.sctenifoldnet already
% does upstream via A0sym/A1sym, and what the sibling ten.i_laplacianembed
% does. On symmetric input it is a bit-exact no-op (verified: max|L_old-L_new|
% is exactly 0), so scTenifoldNet and the GUI callbacks are unaffected; only
% ten.i_knk reaches the changed path.
%
% The off-diagonal blocks are untouched by this: W12 is the identity, hence
% already symmetric, so only the two diagonal blocks move. mu is therefore
% still correct computed above, since sum(0.5*(W1+W1')) == sum(W1).
%
% D keeps the column-sum spelling of the reference (R's
% diag(W) <- -apply(W,2,sum)); on the symmetrized W row and column sums are
% the same thing.
%
% THE SYMMETRIC PATH NEVER FORMS L (2026-10-05). W is two dense gene blocks
% coupled by mu*W12, an identity, so L*x is applied block by block: half
% the multiplications of the dense 2n x 2n product, and no 2n x 2n matrix
% to build. With ndim+1 modes requested rather than 2*ndim (one for the
% null space; 2*ndim is tried next, then the dense fallback), scTenifoldNet's
% alignment at 3000 genes went from 6.4 s to 1.3 s, with drdist unchanged
% to 1e-13 and the gene order identical. The asymmetric path (SYMMETRIZE =
% false, the R reference for Knk) still builds L and calls eigs on it as
% before.
p1 = size(W1, 1);
p2 = size(W2, 1);
if symmetrize
    Lfun = i_operator(W1, W2, mu, laplacian);
    V = i_lomodes(Lfun, p1 + p2, ndim, ...
        @() i_denseL(W1, W2, mu, laplacian, true));
else
    L = i_denseL(W1, W2, mu, laplacian, false);
    V = i_lomodes(L, p1 + p2, ndim, @() L);
end

aln0 = V(1:p1, :);
aln1 = V(p1+1:end, :);
end


function Lfun = i_operator(W1, W2, mu, laplacian)
% L*x for the symmetrized W = [W1s, mu*W12; mu*W12', W2s], without forming
% W or L. W12 = eye(p1, p2) couples gene i of one network with gene i of the
% other, for i up to m = min(p1, p2).
W1 = 0.5 * (W1 + W1.');
W2 = 0.5 * (W2 + W2.');
m = min(size(W1, 1), size(W2, 1));
% Column sums of |W|: the coupling adds |mu| to the first m of each block
d1 = sum(abs(W1), 1).';
d2 = sum(abs(W2), 1).';
d1(1:m) = d1(1:m) + abs(mu);
d2(1:m) = d2(1:m) + abs(mu);

% A KNOCKOUT CHANGES ONE GENE. Under scTenifoldKnk the second network is the
% first with one gene's row zeroed, so after symmetrizing the two blocks
% differ in that gene's row and column only. Each product then reads the
% one dense block once for both halves of x (W1*[x1, x2], one BLAS call)
% and adds the sparse difference, instead of reading two dense blocks:
% ten.knknullcontrol's 60 knockouts at 3000 genes went from 32 s to 23 s.
% Two unrelated networks (scTenifoldNet) differ everywhere and keep the
% two-block product.
W2minusW1 = [];
if isequal(size(W1), size(W2))
    delta = W2 - W1;
    if nnz(delta) <= 0.01*numel(delta)
        W2minusW1 = sparse(delta);
        W2 = [];
    end
    clear delta
end
switch laplacian
    case "normalized"
        s = 1 ./ sqrt([d1; d2]);
        Lfun = @(x) x - s .* i_applyW(W1, W2, W2minusW1, mu, m, s .* x);
    otherwise   % "unnormalized"; checked by the caller
        D = [d1; d2];
        Lfun = @(x) D .* x - i_applyW(W1, W2, W2minusW1, mu, m, x);
end
end


function [aln0, aln1] = i_knockouts(A0, ks, ndim, laplacian)
% One alignment per knocked-out row of A0, sharing A0's preparation. With
% A1 = A0 except row k zeroed, (A1 + 1) differs from (A0 + 1) by -r in row
% k, r = A0(k, :); symmetrized, W2 = W1 - (e_k*r + r'*e_k')/2, so only row
% and column k of the second block, and of its column sums, change.
n = size(A0, 1);
if size(A0, 2) ~= n
    error('ten:i_ma:knockoutsSquare', 'KNOCKOUTS needs a square A0.');
end
W1raw = A0 + 1;
sumW1 = sum(W1raw(:));
W1 = 0.5 * (W1raw + W1raw.');
clear W1raw
d1 = sum(abs(W1), 1).';
nKo = numel(ks);
aln0 = cell(nKo, 1);
aln1 = cell(nKo, 1);
for b = 1:nKo
    k = ks(b);
    r = A0(k, :);
    mu = 0.9 * (2*sumW1 - sum(r)) / (2*n);
    rowk = W1(k, :) - 0.5*r;
    rowk(k) = W1(k, k) - r(k);
    change = rowk - W1(k, :);
    changeCol = change;
    changeCol(k) = 0;   % (k, k) is in the row already
    W2minusW1 = sparse([repmat(k, 1, n), 1:n], [1:n, repmat(k, 1, n)], ...
        [change, changeCol], n, n);
    % Column j ~= k changes only in entry (k, j); column k changes throughout
    d2 = d1 - abs(W1(:, k)) + abs(rowk.');
    d2(k) = sum(abs(rowk));
    d = [d1; d2] + abs(mu);
    switch laplacian
        case "normalized"
            s = 1 ./ sqrt(d);
            Lfun = @(x) x - s .* i_applyW(W1, [], W2minusW1, mu, n, s .* x);
        otherwise   % "unnormalized"; checked by the caller
            Lfun = @(x) d .* x - i_applyW(W1, [], W2minusW1, mu, n, x);
    end
    V = i_lomodes(Lfun, 2*n, ndim, @() i_knockoutDenseL(A0, k, mu, laplacian));
    aln0{b} = V(1:n, :);
    aln1{b} = V(n+1:end, :);
end
end


function L = i_knockoutDenseL(A0, k, mu, laplacian)
% The dense fallback for one knockout
A1 = A0;
A1(k, :) = 0;
L = i_denseL(A0 + 1, A1 + 1, mu, laplacian, true);
end


function y = i_applyW(W1, W2, W2minusW1, mu, m, x)
% W*x for the symmetrized W, block by block. With W2minusW1 given, W2 is
% W1 plus that sparse difference.
p1 = size(W1, 1);
x1 = x(1:p1, :);
x2 = x(p1+1:end, :);
if isempty(W2minusW1)
    y1 = W1 * x1;
    y2 = W2 * x2;
else
    c = size(x1, 2);
    y = W1 * [x1, x2];
    y1 = y(:, 1:c);
    y2 = y(:, c+1:end) + W2minusW1 * x2;
end
y1(1:m, :) = y1(1:m, :) + mu * x2(1:m, :);
y2(1:m, :) = y2(1:m, :) + mu * x1(1:m, :);
y = [y1; y2];
end

function L = i_denseL(W1, W2, mu, laplacian, symmetrize)
% The Laplacian as a matrix, as this function always built it
W12 = eye(size(W1, 2), size(W2, 2));
W = [W1, mu * W12; mu * W12', W2];
if symmetrize
    W = 0.5 * (W + W.');
end
D = sum(abs(W));
switch laplacian
    case "normalized"
        s = 1 ./ sqrt(D);
        L = eye(size(W)) - s.' .* W .* s;
        if symmetrize
            L = 0.5 * (L + L.');    % remove rounding asymmetry so eigs takes Lanczos
        end
    otherwise   % "unnormalized"; checked by the caller
        L = diag(D) - W;
end
end


function V = i_lomodes(L, n, ndim, denseL)
% The ndim eigenvectors of smallest eigenvalue, skipping the null space.
% L is the Laplacian, as a matrix or as a symmetric function handle of size
% n; DENSEL returns it as a matrix, for the fallback.
%
% EIGS IS ITERATIVE AND DOES NOT ALWAYS CONVERGE ON THIS LAPLACIAN. When it
% does not it warns, returns NaN for the eigenvalues that failed, and the
% columns left after the null-space filter can fall below ndim - which used
% to surface as "Index in position 2 exceeds array bounds" from the line
% below, several hours into a scTenifoldNet run. Two things also make a
% partial result unusable rather than merely short: the converged subset is
% not guaranteed to be the smallest ones, and a disconnected aligned graph
% puts more than one eigenvalue at zero.
%
% So: try eigs, asking for one mode more than needed and then, if the null
% space took more, for 2*ndim; on anything less than full convergence redo
% it with a dense symmetric eig, which is exact but grows as n^3, so it is a
% fallback rather than the default. A matrix L keeps the single 2*ndim
% request it always had.
if isa(L, 'function_handle')
    counts = [ndim + 1, 2 * ndim];
else
    counts = 2 * ndim;
end
for k = counts
    [V, d] = i_eigs(L, n, k);
    ok = ~isempty(V) && ~any(isnan(d)) && sum(d >= 1e-8) >= ndim;
    if ok
        break
    end
end
if ~ok
    if isempty(V)
        why = 'eigs failed';
    elseif any(isnan(d))
        why = sprintf('%d of %d eigenvalues did not converge', ...
            sum(isnan(d)), numel(d));
    else
        why = sprintf('only %d eigenvalue(s) above the null-space cutoff', ...
            sum(d >= 1e-8));
    end
    warning('ten:i_ma:denseFallback', ...
        '%s; recomputing the %d smallest modes with a dense eig (n = %d).', ...
        why, ndim, n);
    [V, d] = eig(full(denseL()), 'vector');
    V = real(V);
    [d, ind] = sort(real(d));
    V = V(:, ind);
end

V = V(:, d >= 1e-8);
assert(size(V, 2) >= ndim, 'ten:i_ma:tooFewModes', ...
    ['only %d eigenvalue(s) above 1e-8 but %d modes needed - the aligned ' ...
     'graph is more disconnected than the alignment assumes.'], ...
    size(V, 2), ndim);
V = V(:, 1:ndim);
end


function [V, d] = i_eigs(L, n, k)
% Returns empty V if eigs errors outright, so the caller can fall back on
% that too.
V = []; d = [];
k = min(k, n);
try
    if isa(L, 'function_handle')
        [V, D] = eigs(L, n, k, 'smallestreal', 'IsFunctionSymmetric', true);
    else
        [V, D] = eigs(L, k, 'smallestreal');
    end
catch
    return
end
d = real(diag(D));
[d, ind] = sort(d);
V = real(V(:, ind));
end
