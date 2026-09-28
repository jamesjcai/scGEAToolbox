function [aln0, aln1] = i_ma(A0, A1, ndim, laplacian, symmetrize)
% MA - manifold alignment
%
% [aln0, aln1] = ten.i_ma(A0, A1)
% [aln0, aln1] = ten.i_ma(A0, A1, ndim, laplacian)
% [aln0, aln1] = ten.i_ma(A0, A1, ndim, laplacian, symmetrize)
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

if nargin < 3, ndim = 30; end
if nargin < 4, laplacian = "normalized"; end
if nargin < 5, symmetrize = true; end
mu = 0.9;

W1 = A0 + 1;
W2 = A1 + 1;

W12 = eye(size(W1, 2), size(W2, 2));
mu = mu * (sum(W1(:)) + sum(W2(:))) / (2 * sum(W12(:)));
W = [W1, mu * W12; mu * W12', W2];

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
% D keeps the column-sum spelling rather than the sum(...,2) that
% ten.i_laplacianembed uses. On the now-symmetric W the two are mathematically
% identical, but they accumulate in different orders, and the column form is
% the one the old code used: measured on symmetric input, this branch differs
% from the previous one by exactly 0, where sum(abs(W),2) differs by 1.1e-12.
% That keeps scTenifoldNet's published results bit-identical rather than
% merely equal to rounding.
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
    case "unnormalized"
        L = diag(D) - W;
    otherwise
        error('ten:i_ma:badLaplacian', ...
            'LAPLACIAN must be "normalized" or "unnormalized". Received "%s".', ...
            laplacian);
end

V = i_lomodes(L, ndim);

p1 = size(W1, 1);
aln0 = V(1:p1, :);
aln1 = V(p1+1:end, :);
end


function V = i_lomodes(L, ndim)
% The ndim eigenvectors of smallest eigenvalue, skipping the null space.
%
% EIGS IS ITERATIVE AND DOES NOT ALWAYS CONVERGE ON THIS LAPLACIAN. When it
% does not it warns, returns NaN for the eigenvalues that failed, and the
% columns left after the null-space filter can fall below ndim - which used
% to surface as "Index in position 2 exceeds array bounds" from the line
% below, several hours into a scTenifoldNet run. Two things also make a
% partial result unusable rather than merely short: the converged subset is
% not guaranteed to be the smallest ones, and a disconnected aligned graph
% puts more than ndim eigenvalues at zero, so ndim*2 is not always enough
% cushion.
%
% So: try eigs, and on anything less than full convergence redo it with a
% dense symmetric eig, which is exact and affordable because L is dense
% already. That costs seconds at the usual sizes and grows as n^3, so it is
% a fallback rather than the default.
[V, d] = i_eigs(L, ndim);
if isempty(V) || any(isnan(d)) || sum(d >= 1e-8) < ndim
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
        why, ndim, size(L, 1));
    [V, d] = eig(full(L), 'vector');
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


function [V, d] = i_eigs(L, ndim)
% ndim*2 leaves room for the null space. Returns empty V if eigs errors
% outright, so the caller can fall back on that too.
V = []; d = [];
try
    [V, D] = eigs(L, ndim*2, 'smallestreal');
catch
    return
end
d = real(diag(D));
[d, ind] = sort(d);
V = real(V(:, ind));
end
