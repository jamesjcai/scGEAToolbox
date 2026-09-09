function [T, Zx, Zy] = sc_dpg(X, Y, g, setmatrx, setnames, setgenes, ranknorm, bgsubtract)
% SC_DPG - Differential Program (DP) analysis between two groups of cells
%
% Instead of testing individual genes, this function projects cell expression
% through a gene set membership matrix to obtain per-program activity scores,
% then compares those scores between two groups using a Wilcoxon rank-sum test.
% This is the statistical analogue of gene-set-level differential expression.
%
% Algorithm:
%   1. Intersect the gene set gene list with the measured genes.
%   2. (Optional) Apply within-cell rank normalization or background subtraction.
%   3. Compute program activity scores: Z = setmatrx * X  (sets x cells),
%      where each row of setmatrx is a binary membership vector for one program.
%   4. For each program, apply Wilcoxon rank-sum test between the two groups.
%   5. Compute log2 fold-change of mean activity scores (group1 / group2).
%   6. Adjust p-values with Benjamini-Hochberg FDR correction.
%   7. Return programs with |log2FC| >= 1, adjusted p < 0.01, gene set size >= 5.
%
% Inputs:
%   X, Y        : log-normalized expression matrices (genes x cells) for each group
%   g           : gene name list matching rows of X and Y
%   setmatrx    : binary matrix (programs x genes) — gene set membership
%   setnames    : program name list (length = rows of setmatrx)
%   setgenes    : gene list corresponding to columns of setmatrx
%   ranknorm    : (optional) if true, rank genes within each cell before scoring,
%                 making scores robust to expression scale (slower; default: false).
%                 Inspired by GSVA's rank-based enrichment scoring (PMID:23323831).
%   bgsubtract  : (optional) if true, subtract per-cell mean expression before
%                 scoring to remove global transcriptional activity bias (default: false).
%                 Inspired by GSVA's implicit background correction via within-sample ranking.
%
% Output:
%   T : table of significant differential programs, sorted by adjusted p-value
%
% See also: sc_deg, sc_dvg, pkg.e_getgenesets

% X = log1p(sc_norm(X));
if nargin < 7, ranknorm    = false; end
if nargin < 8, bgsubtract  = false; end
if nargin < 6
    [setmatrx, setnames, setgenes] = pkg.e_getgenesets;
end
[~, ix, iy] = intersect(upper(setgenes), upper(g));

setgenes = setgenes(ix);
setmatrx = setmatrx(:, ix);    % s x g

X = X(iy, :);              % g x c
Y = Y(iy, :);              % g x c
g = g(iy);

% Optional: rank genes within each cell (GSVA-style, slower but scale-robust)
if ranknorm
    X = i_colrank(X);
    Y = i_colrank(Y);
end

% Optional: subtract per-cell mean expression (removes global activity bias)
if bgsubtract
    X = X - mean(X, 1);
    Y = Y - mean(Y, 1);
end

Zx = setmatrx * X;               % s x c
Zy = setmatrx * Y;               % s x c

gsetsize = full(sum(setmatrx ~= 0, 2));   % gene number per set (counts
                               % nonzeros so signed regulons from
                               % pkg.e_getgenesets('TFsigned') are not
                               % reduced to nActivated - nRepressed)

n = size(Zx, 1);

% NAN, not ONES. A program whose genes are all absent from G fails the
% any(setmatrx(k, :)) guard below and is never tested, so its entry keeps
% whatever it was initialised to. A hard-coded 1 is finite, and PKG.E_FDR
% counts every finite entry towards the family size -- its own header says
% "a hypothesis that was never tested is not one of the hypotheses being
% corrected over" -- so these placeholders inflated the adjusted p-value of
% every real program. Measured on a 9-program fixture where 5 programs had
% no measured genes: the surviving program's adjusted p went from 0.008658
% to 0.019481, a factor of 2.25, which moved it across the p_adj < 0.01 cut
% this function applies at the end. The all-zero rows come from the
% unconditional intersection above, so any real gene-set collection
% produces them in bulk.
p_val = nan(n, 1);
avg_log2FC = nan(n, 1);
v1 = nan(n, 1);
v2 = nan(n, 1);
n1 = nan(n, 1);
n2 = nan(n, 1);
m1 = nan(n, 1);
m2 = nan(n, 1);

for k = 1:n
    if any(setmatrx(k, :))
        a = Zx(k, :);
        b = Zy(k, :);
        p_val(k) = ranksum(a, b);

        % The descriptives are computed for every tested program. They used
        % to sit inside `if ~isnan(p_val(k)) && p_val(k) < 1e-3`, and since
        % avg_log2FC starts as NaN and the table assembly below deletes
        % every row with isnan(T.avg_log2FC), that gate silently deleted any
        % program whose raw p-value was above 1e-3 -- including programs
        % that pass all three criteria step 7 of the help lists. It is
        % undocumented: step 5 says the fold change is computed for each
        % program, and nothing mentions a raw-p threshold.
        %
        % It bites hardest where there are few cells, because the rank-sum
        % p-value then has a floor. With six cells per group the smallest
        % two-sided p obtainable is 2/nchoosek(12,6) = 2.164e-3, so NO
        % program could ever clear the gate however complete the
        % separation: measured on a program with that p, |log2FC| = 1.18
        % and 8 genes, sc_dpg returned an empty table.
        %
        % The gate was a vestige guarding the expensive NBINFIT calls just
        % below, which are commented out; what runs now is MEAN, which costs
        % nothing. SC_DEG has no such gate -- it computes the fold change
        % for every gene and applies only the FDR cut.
        %
        % [ax]=nbinfit(a);
        % [bx]=nbinfit(b);
        [ax] = mean(a);
        [bx] = mean(b);
        ratio = ax(1) ./ bx(1);
        if ratio > 0
            % sign(ax-bx) gives correct direction when both values are
            % negative (bgsubtract case): ratio is positive but the larger
            % absolute negative value is actually the lower score.
            avg_log2FC(k) = sign(ax(1)) * log2(ratio);
        end
        % ratio <= 0 means means have opposite signs; avg_log2FC stays NaN
        % and is filtered out downstream — log2 of negative is complex
        v1(k) = ax(1);
        v2(k) = bx(1);
        n1(k) = numel(a);
        n2(k) = numel(b);
        m1(k) = sum(a > 0);
        m2(k) = sum(b > 0);
    end
end
% warning on
% PKG.E_FDR replaces the two-branch block that used to sit here. Its MAFDR
% branch and its PKG.E_FDR_BH branch are the same algorithm on clean input
% -- they agree to 2e-16 -- but they disagree whenever the p-values carry
% NaN, which is what RANKSUM returns for a gene with no counts in either
% group and what a per-cell-type run produces in bulk. One drops those from
% the family, the other counts them, so the same command gave a different
% answer depending on whether the Bioinformatics Toolbox was installed.
p_val_adj = pkg.e_fdr(p_val);

T = table(setnames, gsetsize, v1, v2, avg_log2FC, m1, n1, m2, n2, p_val, p_val_adj);
T(isnan(T.p_val) | isnan(T.avg_log2FC) | abs(T.avg_log2FC) < 1, :) = [];
T = sortrows(T, 'p_val_adj', 'ascend');
T = T(T.p_val_adj < 0.01 & T.gsetsize >= 5, :);
end

function R = i_colrank(A)
% Rank genes within each cell (column). Ties are averaged.
R = zeros(size(A));
for j = 1:size(A, 2)
    R(:, j) = tiedrank(A(:, j));
end
end
