function [score] = sc_cellscore(X, genelist, tgsPos, tgsNeg, methodid)
%SC_CELLSCORE  Cell-level gene signature scoring.
%  score = SC_CELLSCORE(X, genelist, tgsPos, tgsNeg, methodid)
%
%  X         : G x C expression matrix (genes x cells).
%  genelist  : G x 1 string/cell array of gene names.
%  tgsPos    : positive marker genes (string array)
%  tgsNeg    : negative marker genes (string array), can be []. Used by
%              method 2 only; methods 1 and 3 have no negative-marker term
%              and warn if one is supplied.
%  methodid  : 1 = UCell (rank-based, PMID:34285779)
%              2 = AddModuleScore/Seurat (default)
%              3 = AUCell (AUC recovery curve)
%
% see also: PKG.E_CELLSCORES, SC_CELLCYCLESCORE

if nargin < 5 || isempty(methodid), methodid = 2; end
if nargin < 4, tgsNeg = []; end
if nargin < 3 || isempty(tgsPos)
    error('USAGE: >>[score]=sc_cellscore(X,genelist,tgsPos);');
end

if ~any(matches(genelist, tgsPos, 'IgnoreCase', true))
    score = NaN(size(X, 2), 1);
    warning('No feature genes found in GENELIST. NaN scores returned');
    return;
end

% TGSNEG is documented as a general input but only method 2 has anywhere
% to put it: i_ucell and i_aucell take no such argument, so it used to be
% accepted and dropped without a word. Callers that pass a signature's
% negative markers -- pkg.e_cellscores does, from the NegativeMarkers column
% of the shipped table -- got a positive-only score back while believing
% otherwise, and the method is chosen in a dialog well away from the call.
if ~isempty(tgsNeg) && (methodid == 1 || methodid == 3)
    warning('sc_cellscore:negativeMarkersIgnored', ...
        ['Method %d has no negative-marker term, so the %d negative ', ...
        'gene(s) supplied are ignored. Use methodid 2 to subtract ', ...
        'them.'], methodid, numel(string(tgsNeg)));
end

switch methodid
    case 1
        score = i_ucell(X, genelist, tgsPos);
    case 2
        score = i_admdl(X, genelist, tgsPos, tgsNeg);
    case 3
        score = i_aucell(X, genelist, tgsPos);
    otherwise
        error('Unknown methodid %d. Use 1 (UCell), 2 (AddModuleScore), or 3 (AUCell).', methodid);
end

end


%% ---- Method 1: UCell (rank-based) ----
function [score] = i_ucell(X, genelist, tgsPos, maxRank)
% UCell rank-based signature scoring (Andreatta & Carmona, 2021).
% ref: https://doi.org/10.1016/j.csbj.2021.06.043
% ref: https://github.com/carmonalab/UCell

if nargin < 4 || isempty(maxRank), maxRank = 1500; end

idx = matches(genelist, tgsPos, 'IgnoreCase', true);
n = sum(idx);

% Per-cell ranks with highest expression at rank 1; average ties.
R = tiedrank(-full(X));
R(R > maxRank) = maxRank + 1;

% Mann-Whitney U statistic per cell: rank sum minus its minimum n(n+1)/2.
rankSum = sum(R(idx, :), 1);
u = rankSum - (n * (n + 1)) / 2;
score = 1 - u / (n * maxRank);

% Cells whose signature genes all rank beyond maxRank score 0 (UCell).
score(all(R(idx, :) > maxRank, 1)) = 0;
score = score(:);
end


%% ---- Method 2: AddModuleScore/Seurat ----
function [score] = i_admdl(X, genelist, tgsPos, tgsNeg, nbin, ctrl)
% AddModuleScore - Seurat-style scoring
% ref: https://github.com/satijalab/seurat/blob/master/R/utilities.R
% ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8271111/

if nargin < 6, ctrl = 100; end
if nargin < 5, nbin = 24; end
if nargin < 4, tgsNeg = []; end

if issparse(X), X = full(X); end

X = sc_norm(X);
X = log1p(X);

[score] = i_admdl_calculate(X, genelist, tgsPos, 1, nbin, ctrl);
if ~isempty(tgsNeg) && any(strlength(tgsNeg) > 0)
    [s] = i_admdl_calculate(X, genelist, tgsNeg, -1, nbin, ctrl);
    score = score + s;
end
end

function [score] = i_admdl_calculate(X, genelist, tgs, directtag, nbin, ctrl)
if nargin < 6, ctrl = 100; end
if nargin < 5, nbin = 24; end
if nargin < 4, directtag = 1; end

cluster_length = size(X, 1);
data_avg = mean(X, 2);

% Break exact ties (e.g. all-zero genes) before binning, matching Seurat's
% addition of rnorm(n)/1e30 to data.avg prior to cut_number.
data_avg = data_avg + randn(size(data_avg)) / 1e30;

[~, I] = sort(data_avg);
data_avg = data_avg(I);
gsorted = genelist(I);
Xsorted = X(I, :);

% Equal-frequency bins over sorted mean expression (Seurat cut_number).
assigned_bin = zeros(cluster_length, 1);
bin_size = cluster_length / nbin;
for i = 1:nbin
    bin_match = data_avg <= data_avg(round(bin_size*i));
    pos_avail = (assigned_bin == 0);
    assigned_bin(pos_avail & bin_match) = i;
end

idx = matches(gsorted, tgs, 'IgnoreCase', true);

% Draw ctrl control genes from each feature gene's own bin (per-feature
% sampling, matching Seurat) rather than from a pool of all feature bins.
ctrl_cell = cell(length(tgs), 1);
for i = 1:length(tgs)
    gi = find(matches(gsorted, tgs(i), 'IgnoreCase', true), 1);
    if isempty(gi)
        continue;
    end
    bin_genes = gsorted(assigned_bin == assigned_bin(gi));
    k = min(ctrl, numel(bin_genes));
    ctrl_cell{i} = bin_genes(randsample(numel(bin_genes), k));
end
ctrl_use = unique(vertcat(ctrl_cell{:}));

ctrl_score = mean(Xsorted(matches(gsorted, ctrl_use, 'IgnoreCase', true), :), 1);
features_score = mean(Xsorted(idx, :), 1);

if directtag > 0
    score = transpose(features_score-ctrl_score);
else
    score = transpose(ctrl_score-features_score);
end
end


%% ---- Method 3: AUCell (AUC recovery curve) ----
function [score] = i_aucell(X, genelist, tgsPos, aucMaxRank)
% AUCell - Area Under the recovery Curve scoring (Aibar et al., 2017).
% ref: https://doi.org/10.1038/nmeth.4463
% ref: https://bioconductor.org/packages/AUCell

[nGenes, nCells] = size(X);

% Default aucMaxRank: top 5% of genes (Bioconductor AUCell default).
if nargin < 4 || isempty(aucMaxRank)
    aucMaxRank = ceil(0.05 * nGenes);
end

idx = matches(genelist, tgsPos, 'IgnoreCase', true);
nSet = sum(idx);
if nSet == 0
    score = NaN(nCells, 1);
    return;
end

if issparse(X), X = full(X); end

% Maximum attainable area. At most AUCMAXRANK of the set's genes can sit
% inside the recovery window, so the best case is MIN(nSet, aucMaxRank) of
% them occupying ranks 1, 2, ... -- which is what the MIN is for.
%
% This was nSet*aucMaxRank - nSet*(nSet + 1)/2, correct only while
% nSet <= aucMaxRank. Beyond that it understates the maximum, reaches
% exactly zero at nSet = 2*aucMaxRank - 1, and goes negative after it. On a
% 2000-gene matrix (aucMaxRank = 100) with a signature planted at the very
% top of every cell: nSet = 50 scored 1.000 as it should, nSet = 150 scored
% +1.347 -- an area under a curve cannot exceed 1 -- and nSet = 250 scored
% -0.776, so the strongest signal obtainable came back strongly negative.
% 35 of the 145 signatures shipped in pkg.e_cellscores are larger than 100
% genes and the largest is 1313, so this is an ordinary case, not a corner.
nInWindow = min(nSet, aucMaxRank);
maxArea = nInWindow * aucMaxRank - nInWindow * (nInWindow + 1) / 2;
if maxArea <= 0
    error('sc_cellscore:degenerateRecoveryWindow', ...
        ['A recovery window of %d rank(s) leaves no area to normalise ', ...
        'against. Raise aucMaxRank.'], aucMaxRank);
end

% TIEDRANK, not sort-order. MATLAB's sort is stable, so
%   [~, order] = sort(X(:, c), 'descend'); ranks(order) = 1:nGenes;
% gave every undetected gene a distinct rank in ascending ROW order, and
% the recovery window then admitted whichever zeros happened to sit early
% in GENELIST. Measured on 600 genes with only 5 detected and
% aucMaxRank = 30: a 20-gene signature drawn from rows 1-20 scored 0.7436
% while one drawn from rows 200-219 scored 0.0000, both being entirely
% undetected in every cell. Averaging ties gives every tied gene the same
% rank, so the score no longer depends on gene order; i_ucell above already
% does exactly this.
ranks = tiedrank(-X);              % nGenes x nCells, highest expression first

setRanks = ranks(idx, :);          % nSet x nCells
inWindow = setRanks <= aucMaxRank;

% Area under the step recovery curve equals sum(aucMaxRank - rank_i) over
% the set genes inside the window.
area = sum(inWindow, 1)*aucMaxRank - sum(setRanks.*inWindow, 1);
score = (area./maxArea)';
end
