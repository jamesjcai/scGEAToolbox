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

% Per-cell ranks with highest expression at rank 1, ties averaged -- the
% signature's rows of tiedrank(-X), without ranking the whole matrix. See
% I_SETRANKS.
R = i_setranks(X, idx);
R(R > maxRank) = maxRank + 1;

% Mann-Whitney U statistic per cell: rank sum minus its minimum n(n+1)/2.
rankSum = sum(R, 1);
u = rankSum - (n * (n + 1)) / 2;
score = 1 - u / (n * maxRank);

% Cells whose signature genes all rank beyond maxRank score 0 (UCell).
score(all(R > maxRank, 1)) = 0;
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

% Kept sparse: SC_NORM and LOG1P both preserve sparsity, and only the rows
% of the signature and its control genes are ever averaged. This used to
% FULL() the matrix and then build a gene-sorted dense copy of it as well,
% twice for a signature with negative markers.
X = sc_norm(X);
X = log1p(X);

[score] = i_admdl_calculate(X, genelist, tgsPos, 1, nbin, ctrl);
% Apply the negative term only when at least one of its genes is present in
% the data. A tgsNeg whose genes are all absent is an ordinary situation - a
% masking enzyme simply not expressed in this tissue - and must leave the
% positive score untouched. Without this guard the control-gene pool comes
% back empty and MATCHES errors on it.
if ~isempty(tgsNeg) && any(strlength(string(tgsNeg)) > 0) && ...
        any(matches(genelist, tgsNeg, 'IgnoreCase', true))
    [s] = i_admdl_calculate(X, genelist, tgsNeg, -1, nbin, ctrl);
    score = score + s;
end
end

function [score] = i_admdl_calculate(X, genelist, tgs, directtag, nbin, ctrl)
if nargin < 6, ctrl = 100; end
if nargin < 5, nbin = 24; end
if nargin < 4, directtag = 1; end

cluster_length = size(X, 1);
data_avg = full(mean(X, 2));

% Break exact ties (e.g. all-zero genes) before binning, matching Seurat's
% addition of rnorm(n)/1e30 to data.avg prior to cut_number.
data_avg = data_avg + randn(size(data_avg)) / 1e30;

[~, I] = sort(data_avg);
data_avg = data_avg(I);
gsorted = genelist(I);
% No XSORTED = X(I, :): a permuted copy of the whole matrix, when only a few
% hundred rows are read. Rows are picked below as I(mask), which is the
% same rows in the same order.

% Equal-frequency bins over sorted mean expression (Seurat cut_number).
%
% MAX(1, ...) because ROUND(bin_size*i) is 0 for the first bins whenever
% there are fewer genes than bins: with 10 genes and nbin 24 the bin size is
% 0.42 and the very first iteration indexed data_avg(0). That is a legal
% input - a small panel, or a matrix subset to a handful of genes - and it
% crashed with "Array indices must be positive integers" from three frames
% down. For any matrix with at least nbin genes ROUND(bin_size*i) is already
% 1 or more, so this clamp changes nothing for existing callers; below that
% the leading bins collapse onto the lowest-expressed gene, which is the
% sensible degenerate behaviour.
assigned_bin = zeros(cluster_length, 1);
bin_size = cluster_length / nbin;
for i = 1:nbin
    bin_match = data_avg <= data_avg(max(1, round(bin_size*i)));
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

% None of TGS is in the data: there is no signal and no matched background,
% so the score is undefined rather than zero. VERTCAT of empty cells gives a
% 0-by-0 double, which MATCHES would reject, so bail before it.
if isempty(ctrl_use)
    score = NaN(size(X, 2), 1);
    return;
end

ctrl_score = full(mean(X(I(matches(gsorted, ctrl_use, 'IgnoreCase', true)), :), 1));
features_score = full(mean(X(I(idx), :), 1));

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
% I_SETRANKS returns exactly the set's rows of tiedrank(-X).
setRanks = i_setranks(X, idx);     % nSet x nCells, highest expression first
inWindow = setRanks <= aucMaxRank;

% Area under the step recovery curve equals sum(aucMaxRank - rank_i) over
% the set genes inside the window.
area = sum(inWindow, 1)*aucMaxRank - sum(setRanks.*inWindow, 1);
score = (area./maxArea)';
end


%% ---- per-cell ranks of the signature genes ----
function R = i_setranks(X, idx)
% R = the rows IDX of tiedrank(-full(X)): each cell's genes ranked from the
% highest value down, ties averaged. Only the signature's rows are
% returned, so R is nSet x nCells.
%
% Ranking the whole matrix, as UCell and AUCell did, densified it and sorted
% every gene of every cell: 11 s and ~20 GB at 20000 x 50000 cells. Here each
% cell's stored nonzeros are sorted once; its zeros are one tied block,
% ranked after the positives, so a signature gene that is zero in a cell
% takes that block's average rank. A value that is NaN anywhere sends the
% matrix to TIEDRANK itself, whose NaN handling this does not reproduce.
if any(isnan(nonzeros(X)))
    R = tiedrank(-full(X));
    R = R(idx, :);
    return;
end

[G, C] = size(X);
setpos = zeros(G, 1);
setpos(idx) = 1:nnz(idx);                  % gene row -> row of R
R = zeros(nnz(idx), C);

perCol = max(1, nnz(X)/max(C, 1));
step = max(1, floor(1e7/perCol));
for c0 = 1:step:C
    cols = c0:min(C, c0 + step - 1);
    [i, j, v] = find(X(:, cols));
    i = i(:);
    j = j(:);
    v = full(double(v(:)));

    % Each cell's nonzeros from the highest down.
    [~, o] = sortrows([j, -v]);
    i = i(o);
    j = j(o);
    v = v(o);

    nc = numel(cols);
    nnzc = accumarray(j, 1, [nc, 1]);
    npos = accumarray(j, v > 0, [nc, 1]);
    n0 = G - nnzc;

    % Position among the cell's nonzeros, averaged over ties.
    first = cumsum([0; nnzc(1:end-1)]);
    pos = (1:numel(i)).' - first(j);
    newtie = [true; diff(j) ~= 0 | diff(v) ~= 0];
    tid = cumsum(newtie);
    tsize = accumarray(tid, 1);
    tfirst = pos(newtie);
    avgpos = tfirst(tid) + (tsize(tid) - 1)/2;

    % Positives rank first, then the zero block, then any negatives.
    r = avgpos + (v < 0).*n0(j);
    zeroRank = npos + (n0 + 1)/2;

    Rb = repmat(zeroRank.', size(R, 1), 1);
    keep = setpos(i) > 0;
    Rb(sub2ind(size(Rb), setpos(i(keep)), j(keep))) = r(keep);
    R(:, cols) = Rb;
end
end
