function [T] = sc_dpg(X, Y, g, setmatrx, setnames, setgenes, ranknorm, bgsubtract)
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

gsetsize = sum(setmatrx, 2);   % gene number per set

n = size(Zx, 1);

p_val = ones(n, 1);
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
        if ~isnan(p_val(k)) && p_val(k) < 1e-3
            % [ax]=nbinfit(a);
            % [bx]=nbinfit(b);
            [ax] = mean(a);
            [bx] = mean(b);
            avg_log2FC(k) = log2(ax(1) ./ bx(1));
            v1(k) = ax(1);
            v2(k) = bx(1);
            n1(k) = numel(a);
            n2(k) = numel(b);
            m1(k) = sum(a > 0);
            m2(k) = sum(b > 0);
        end
    end
end
% warning on
if exist('mafdr.m', 'file')
    p_val_adj = mafdr(p_val, 'BHFDR', true);
else
    [~, ~, ~, p_val_adj] = pkg.e_fdr_bh(p_val);
end

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
