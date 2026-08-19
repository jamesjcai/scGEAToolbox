function [score] = sc_cellscore(X, genelist, tgsPos, tgsNeg, methodid)
%SC_CELLSCORE  Cell-level gene signature scoring.
%  score = SC_CELLSCORE(X, genelist, tgsPos, tgsNeg, methodid)
%
%  X         : G x C expression matrix (genes x cells).
%  genelist  : G x 1 string/cell array of gene names.
%  tgsPos    : positive marker genes (string array)
%  tgsNeg    : negative marker genes (string array), can be []
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

% Maximum attainable area: whole gene set occupies the top ranks.
maxArea = nSet * aucMaxRank - nSet * (nSet + 1) / 2;

score = zeros(nCells, 1);
for c = 1:nCells
    [~, order] = sort(X(:, c), 'descend');
    ranks = zeros(nGenes, 1);
    ranks(order) = 1:nGenes;

    % Ranks of the gene set that fall within the recovery window.
    setRanks = ranks(idx);
    setRanks = setRanks(setRanks <= aucMaxRank);

    % Area under the step recovery curve equals sum(aucMaxRank - rank_i).
    area = numel(setRanks) * aucMaxRank - sum(setRanks);
    score(c) = area / maxArea;
end
end
