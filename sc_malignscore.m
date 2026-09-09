function [score, label, info] = sc_malignscore(cnv, refidx, opts)
%SC_MALIGNSCORE  Score and label cells as malignant from a CNV profile.
%
%   [score, label] = SC_MALIGNSCORE(cnv, refidx) turns the per-gene copy
%   number profile from SC_INFERCNV into one number per cell and, when the
%   distribution of that number is bimodal, splits the cells into malignant
%   and non-malignant.
%
%   The score is the mean squared deviation of a cell's CNV values from 1,
%   i.e. how far its genome departs from the reference karyotype, no matter
%   in which direction or at which locus. A diploid cell scores near zero;
%   a cell carrying several arm-level events scores high.
%
%   USAGE:
%     isref = sce.c_cell_type_tx == "T cells";
%     cnv = sc_infercnv(sce.X, sce.g, isref);
%     [score, label] = sc_malignscore(cnv, isref);
%
%     % sharper separation, using the cell-cell graph to pool neighbours
%     A = sc_knngraph(sce.s, 10);
%     [score, label] = sc_malignscore(cnv, isref, Adjacency=A);
%
%   INPUTS:
%     cnv    - nGenes-by-nCells relative copy number from SC_INFERCNV.
%     refidx - reference cells, as a logical mask over columns or a vector
%              of indices. Use the same set given to SC_INFERCNV: the
%              reference scores calibrate the scale and anchor the split.
%
%   NAME-VALUE:
%     Adjacency - nCells-by-nCells cell-cell similarity or adjacency, e.g.
%                 from SC_KNNGRAPH. When supplied, each cell's CNV profile
%                 is replaced by a weighted average of itself (weight 0.5)
%                 and its neighbours (0.5 shared between them) before
%                 scoring. Dropout makes single-cell CNV profiles noisy,
%                 and cells of one clone share their karyotype, so pooling
%                 neighbours recovers signal that per-cell scoring loses.
%     Normalize - rescale scores so the reference cells span roughly 0 to 1
%                 (default true). Uses the reference 0.5th and 99.5th
%                 percentiles, so a malignant cell lands well above 1.
%     Verbose   - print a short summary (default true).
%
%   OUTPUTS:
%     score - nCells-by-1 malignancy score.
%     label - nCells-by-1 categorical, "malignant" or "nonMalignant". When
%             the scores are unimodal no threshold exists and every cell is
%             labelled nonMalignant; read INFO.HASSPLIT to tell that apart
%             from a genuine absence of malignant cells.
%     info  - struct with the threshold applied, both candidate thresholds
%             it was chosen from (thresholdAll, thresholdObs), whether a
%             split was supported, the reference scores, and the densities
%             behind the decision.
%
%   NOTE. A returned LABEL is not by itself proof of malignancy. Read it
%   together with INFO.HasSplit and the score histogram: when HasSplit is
%   false every cell got the same label by default, not by evidence. Cells
%   with high scores and no epithelial identity are usually a reference
%   problem rather than a discovery.
%
% REF: Guo et al. (2021) Brief Bioinform 22:bbaa127 (scCancer).
%
% See also: SC_INFERCNV, PKG.E_BIMODALTHRES, SC_KNNGRAPH.

arguments
    cnv {mustBeNumeric, mustBeNonempty}
    refidx {mustBeNonempty}
    opts.Adjacency = []
    opts.Normalize (1, 1) logical = true
    opts.Verbose (1, 1) logical = true
end

nCells = size(cnv, 2);
isref = i_torefmask(refidx, nCells);
if ~any(isref)
    error('No reference cells were selected by REFIDX.');
end
if all(isref)
    error(['Every cell is marked as reference, leaving nothing to score. ' ...
        'REFIDX must mark only the normal cells.']);
end

% Pool each cell with its neighbours before scoring, if a graph is given.
if ~isempty(opts.Adjacency)
    W = i_poolweights(opts.Adjacency, nCells);
    scoreAll = i_score(cnv * W);
else
    scoreAll = i_score(cnv);
end

refscore = scoreAll(isref);
if opts.Normalize
    lo = quantile(refscore, 0.005);
    hi = quantile(refscore, 0.995);
    if hi > lo
        scoreAll = (scoreAll - lo) / (hi - lo);
        refscore = scoreAll(isref);
    end
end
score = scoreAll;

% Threshold over every cell, reference included. The reference cells are
% what anchor the non-malignant mode, so a valley found with them present
% is the better-determined one.
%
% scCancer instead cuts at a threshold derived from the observation cells
% alone, because there "Reference" is a separate bundled dataset and
% "Observation" is the whole sample, reliably containing both populations.
% Here the reference is a subset of the sample, so the observation cells
% can be entirely malignant -- select the T cells as reference in a
% tumour-rich sample and that is exactly what happens. Any spurious split
% inside that single cloud is then used to relabel most of the tumour as
% normal. Measured on such a sample: the score separated perfectly
% (AUC 1.000) while the scCancer rule scored 0.437.
obs = scoreAll(~isref);
[allThres, allInfo] = pkg.e_bimodalthres(scoreAll);
[obsThres, obsInfo] = pkg.e_bimodalthres(obs);

thres = allThres;
if isempty(thres)
    thres = obsThres;   % no reference-anchored valley; fall back
end
hasSplit = ~isempty(thres);

lab = repmat("nonMalignant", nCells, 1);
if hasSplit
    lab(scoreAll >= thres) = "malignant";
end
lab(isref) = "nonMalignant";
label = categorical(lab, {'nonMalignant', 'malignant'});

info = struct('threshold', thres, 'thresholdAll', allThres, ...
    'thresholdObs', obsThres, 'hasSplit', hasSplit, 'refScore', refscore, ...
    'isRef', isref(:), 'density', obsInfo, 'densityAll', allInfo, ...
    'pooled', ~isempty(opts.Adjacency));

if opts.Verbose
    if hasSplit
        fprintf(['sc_malignscore: %d of %d non-reference cells called ' ...
            'malignant (threshold %.3f).\n'], ...
            sum(label(~isref) == "malignant"), sum(~isref), thres);
    else
        fprintf(['sc_malignscore: scores are unimodal; no split supported, ' ...
            'all cells labelled nonMalignant.\n']);
    end
end
end


function s = i_score(cnv)
% Mean squared departure from copy-number neutral.
s = sum((double(cnv) - 1).^2, 1).' / size(cnv, 1);
end


function W = i_poolweights(A, nCells)
% Column-stochastic pooling matrix: 0.5 on self, 0.5 spread over neighbours,
% or 1.0 on self for a cell with no neighbours.
if ~isequal(size(A), [nCells, nCells])
    error('Adjacency must be %d-by-%d to match the columns of CNV.', ...
        nCells, nCells);
end

if islogical(A)
    N = A;
else
    % A weighted similarity graph is thresholded to keep roughly ten
    % neighbours per cell, as scCancer does with the Seurat SNN graph.
    nz = nonzeros(A);
    if isempty(nz)
        N = false(nCells);
    elseif isscalar(unique(nz))
        N = A ~= 0;
    else
        N = A > quantile(nz, max(0, 1 - nCells * 10 / numel(nz)));
    end
end
N = double(N);
N(1:nCells+1:end) = 0;   % a cell is not its own neighbour

nNeighbor = sum(N, 2);
share = zeros(nCells, 1);
share(nNeighbor > 0) = 0.5 ./ nNeighbor(nNeighbor > 0);
selfWeight = repmat(0.5, nCells, 1);
selfWeight(nNeighbor == 0) = 1;

% SPDIAGS, not DIAG: a kNN adjacency is sparse and stays sparse through the
% scaling, but DIAG(selfWeight) would materialise a full nCells-by-nCells
% matrix -- 3.2 GB at 20k cells -- and densify the sum.
W = (N .* share + spdiags(selfWeight, 0, nCells, nCells)).';
end


function isref = i_torefmask(refidx, nCells)
if islogical(refidx)
    if numel(refidx) ~= nCells
        error('REFIDX is a logical mask of length %d but CNV has %d columns.', ...
            numel(refidx), nCells);
    end
    isref = refidx(:).';
else
    if any(refidx < 1 | refidx > nCells | mod(refidx, 1) ~= 0)
        error('REFIDX holds indices outside 1..%d.', nCells);
    end
    isref = false(1, nCells);
    isref(refidx) = true;
end
end
