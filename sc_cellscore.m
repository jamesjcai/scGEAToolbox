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
function [score] = i_ucell(X, genelist, tgsPos)
% UCell-inspired rank-based scoring
% https://doi.org/10.1016/j.csbj.2021.06.043

genelist = upper(genelist);
tgsPos = upper(tgsPos);

idx1 = matches(genelist, tgsPos, 'IgnoreCase', true);
n1 = sum(idx1);

if issparse(X)
    try
        X = full(X);
    catch
        warning('Could not convert sparse matrix to full. Proceeding with sparse.');
    end
end

R = tiedrank(-X);
R(R > 1500) = 1500 + 1;
u = sum(R(idx1, :)) - (n1 * (n1 - 1)) / 2;
score = 1 - u / (n1 * 1500);
score = score(:);
end


%% ---- Method 2: AddModuleScore/Seurat ----
function [score] = i_admdl(X, genelist, tgsPos, tgsNeg, nbin, ctrl)
% AddModuleScore - Seurat-style scoring
% ref: https://github.com/satijalab/seurat/blob/master/R/utilities.R
% ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8271111/

if nargin < 6, ctrl = 5; end
if nargin < 5, nbin = 25; end
if nargin < 4, tgsNeg = []; end

try
    if issparse(X), X = full(X); end
catch
    warning('Keep using sparse X.');
end

X = sc_norm(X);
X = log1p(X);

[score] = i_admdl_calculate(X, genelist, tgsPos, 1, nbin, ctrl);
if ~isempty(tgsNeg) && any(strlength(tgsNeg) > 0)
    [s] = i_admdl_calculate(X, genelist, tgsNeg, -1, nbin, ctrl);
    score = score + s;
end
end

function [score] = i_admdl_calculate(X, genelist, tgs, directtag, nbin, ctrl)
if nargin < 6, ctrl = 5; end
if nargin < 5, nbin = 25; end
if nargin < 4, directtag = 1; end

cluster_length = size(X, 1);
data_avg = mean(X, 2);
[~, I] = sort(data_avg);

data_avg = data_avg(I);
gsorted = genelist(I);
Xsorted = X(I, :);

assigned_bin = zeros(cluster_length, 1);
bin_size = cluster_length / nbin;
for i = 1:nbin
    bin_match = data_avg <= data_avg(round(bin_size*i));
    pos_avail = (assigned_bin == 0);
    assigned_bin(pos_avail & bin_match) = i;
end

idx = matches(gsorted, tgs, 'IgnoreCase', true);
selected_bins = unique(assigned_bin(idx));
samebin_genes = gsorted(ismember(assigned_bin, selected_bins));
ctrl_use = [];
for i = 1:length(tgs)
    ctrl_use = [ctrl_use; ...
        randsample(samebin_genes, ctrl)];
end
ctrl_use = unique(ctrl_use);

ctrl_score = mean(Xsorted(matches(gsorted, ctrl_use, 'IgnoreCase', true), :), 1);
features_score = mean(Xsorted(idx, :), 1);

if directtag > 0
    score = transpose(features_score-ctrl_score);
else
    score = transpose(ctrl_score-features_score);
end
end


%% ---- Method 3: AUCell (AUC recovery curve) ----
function [score] = i_aucell(X, genelist, tgsPos)
% AUCell - Area Under the recovery Curve scoring

[nGenes, nCells] = size(X);

% Default aucMaxRank: 5% of genes, minimum 50
aucMaxRank = max(50, round(0.05 * nGenes));

% Find gene indices for the positive markers
idx = matches(genelist, tgsPos, 'IgnoreCase', true);
geneIndices = find(idx);

if isempty(geneIndices)
    score = NaN(nCells, 1);
    return;
end

% Create rankings for each cell (rank by descending expression)
rankings = zeros(nGenes, nCells);
for i = 1:nCells
    [~, rankings(:,i)] = sort(X(:,i), 'descend');
end

% Calculate AUC for each cell
score = zeros(nCells, 1);
for cellIdx = 1:nCells
    score(cellIdx) = i_aucell_auc(rankings(:, cellIdx), geneIndices, aucMaxRank);
end
end

function auc = i_aucell_auc(geneRanks, geneSet, maxRank)
% Calculate AUC for a single gene set in a single cell

setRanks = geneRanks(geneSet);
setRanks = setRanks(setRanks <= maxRank);

if isempty(setRanks)
    auc = 0;
    return;
end

setRanks = sort(setRanks);
nGenesInSet = length(setRanks);

x = [0; setRanks; maxRank];
y = [0; (1:nGenesInSet)'/nGenesInSet; 1];

auc = trapz(x, y) / maxRank;
auc = max(0, (auc - 0.5) * 2);
end
