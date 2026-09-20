function [cnv, T, seg] = sc_infercnv(X, genelist, refidx, opts)
%SC_INFERCNV  Infer large-scale copy number variation from scRNA-seq.
%
%   cnv = SC_INFERCNV(X, genelist, refidx) estimates a relative copy number
%   profile for every cell by taking a moving average of expression along
%   each chromosome and contrasting it with a set of reference cells that
%   are assumed to be karyotypically normal.
%
%   The premise is that a chromosomal gain or loss shifts the expression of
%   every gene in the affected region by a small, consistent amount. No
%   single gene carries that signal above its own biological and technical
%   noise, but a running mean over ~100 neighbouring genes does, because
%   the gene-specific noise averages out while the regional shift does not.
%
%   USAGE:
%     isref = sce.c_cell_type_tx == "T cells";
%     [cnv, T] = sc_infercnv(sce.X, sce.g, isref);
%     [score, label] = sc_malignscore(cnv, isref);
%
%   INPUTS:
%     X        - nGenes-by-nCells matrix of raw UMI counts, full or sparse.
%     genelist - nGenes-by-1 gene symbols, or Ensembl gene IDs (detected
%                automatically).
%     refidx   - reference cells, as a logical mask over the columns of X
%                or a vector of column indices. These should be normal
%                cells of the same sample, most often T/NK or myeloid
%                cells. Pick a type that is abundant and unlikely to be
%                aneuploid; the whole result is expressed relative to it.
%
%   NAME-VALUE:
%     Genome       - 'hg38' (default), 'hg19' or 'mm10'. Gene coordinates
%                    come from PKG.E_GENECOORDINATES.
%     MinExprMean  - drop genes whose mean count over all cells is below
%                    this (default 0.1). Low-expression genes contribute
%                    noise but no regional signal.
%     MinCells     - drop genes detected in fewer than this many cells
%                    (default 3).
%     WindowLength - genes per moving average, forced odd (default 101).
%                    Larger windows suppress noise but blur breakpoints and
%                    can erase focal events.
%     SdAmplifier  - half-width, in reference standard deviations, of the
%                    dead zone flattened to no-change (default 1.0). Larger
%                    values report fewer but more confident alterations;
%                    0 disables the denoising step.
%     Segment      - also fit piecewise-constant segments to the profile
%                    (default false; implied by asking for the third
%                    output). This is the step copykat and SCEVAN add on
%                    top of the moving average, fitted jointly across
%                    cells so that every cell ends up on the same
%                    boundaries and the columns stay comparable.
%
%                    It is for reporting, not for detection. On simulated
%                    focal gains spanning 1.08x to 1.6x over 60 to 150
%                    genes, calling altered genes from the segmented
%                    profile matched calling them from CNV itself in seven
%                    of eight conditions and was slightly worse in the
%                    eighth: the moving average has already imposed the
%                    smoothness that segmentation would otherwise supply.
%                    What it does give is SEG.TABLE -- a few hundred
%                    regions with coordinates, instead of a value per gene
%                    -- which is what you want to export, or to intersect
%                    with a gene of interest. Read a gain off CNV; describe
%                    where it starts and stops from SEG.
%     SegStopRatio - how hard to merge before stopping, as a multiple of
%                    the noise level (default 2). Lower gives more, shorter
%                    segments. Used only when segmenting.
%     UseSingle    - hold the working matrix in single precision
%                    (default true). CNV values are ratios near 1, so
%                    single is ample, and it halves peak memory.
%     Verbose      - print progress (default true).
%
%   OUTPUTS:
%     cnv - nGenesKept-by-nCells relative copy number, centred on 1. Values
%           above 1 suggest a gain, below 1 a loss. Rows are ordered along
%           the genome, so IMAGESC(cnv') draws the familiar heat map with
%           chromosomes running left to right.
%     T   - table of the retained genes in row order, with the variables
%           Gene, EnsemblID, Chr, Start and Stop.
%     seg - segmentation, when asked for: a struct with .cnv, the same size
%           as CNV but constant within each segment, .breakpoints, the rows
%           where segments start, and .table, one row per segment carrying
%           its chromosome, first and last gene, and span in bases. Empty
%           when Segment is false.
%
%   The pipeline follows scCancer's implementation of the inferCNV
%   algorithm: library-size normalisation, Anscombe variance stabilisation,
%   log2 transform, symmetric clipping at the mean per-cell extreme,
%   per-chromosome triangular-kernel moving average, per-cell median
%   centring, subtraction of the reference profile, inversion out of log
%   space, flattening of the reference dead zone, and a final outlier
%   clamp.
%
%   NOTE ON REFERENCE CHOICE. Everything here is relative. If the reference
%   cells are themselves aneuploid, or are so few that their mean profile
%   is noisy, the shared component cancels and real events disappear. Use a
%   few hundred reference cells where possible, and prefer immune cells
%   over stromal ones. Sex chromosomes are excluded because their copy
%   number is confounded by donor sex.
%
% REF: Patel et al. (2014) Science 344:1396 (the moving-average idea);
%      Tickle et al., inferCNV of the Trinity CTAT Project;
%      Guo et al. (2021) Brief Bioinform 22:bbaa127 (scCancer).
%
% See also: SC_MALIGNSCORE, PKG.E_GENECOORDINATES, RUN.R_INFERCNV,
%           RUN.ML_SCEVAN

arguments
    X {mustBeNumeric, mustBeNonempty}
    genelist (:, 1) string
    refidx {mustBeNonempty}
    opts.Genome (1, :) char {mustBeMember(opts.Genome, {'hg38', 'hg19', 'mm10'})} = 'hg38'
    opts.MinExprMean (1, 1) double {mustBeNonnegative} = 0.1
    opts.MinCells (1, 1) double {mustBeNonnegative} = 3
    opts.WindowLength (1, 1) double {mustBePositive} = 101
    opts.SdAmplifier (1, 1) double {mustBeNonnegative} = 1.0
    opts.Segment (1, 1) logical = false
    opts.SegStopRatio (1, 1) double {mustBePositive} = 2
    opts.UseSingle (1, 1) logical = true
    opts.Verbose (1, 1) logical = true
end

[nGenes, nCells] = size(X);
if numel(genelist) ~= nGenes
    error('GENELIST has %d entries but X has %d rows. They must match.', ...
        numel(genelist), nGenes);
end

isref = i_torefmask(refidx, nCells);
if ~any(isref)
    error(['No reference cells were selected. REFIDX must mark at least ' ...
        'one column of X as a karyotypically normal cell.']);
end
if opts.Verbose
    fprintf('sc_infercnv: %d cells, %d of them reference.\n', nCells, sum(isref));
end

% ---- Map genes to genomic coordinates and order them along the genome ----
[tf, loc, T] = pkg.e_matchgenecoords(genelist, opts.Genome);

X = X(tf, :);
T = T(loc(tf), :);
[~, ord] = sortrows([T.Chr, T.Start]);
X = X(ord, :);
T = T(ord, :);
if opts.Verbose
    fprintf('  %d of %d genes placed on chromosomes 1-%d.\n', ...
        height(T), nGenes, max(T.Chr));
end

% ---- Drop genes too sparse to contribute to a regional average ----
keep = full(mean(X, 2)) >= opts.MinExprMean & full(sum(X > 0, 2)) >= opts.MinCells;
if ~any(keep)
    error(['Every gene was filtered out. Lower MinExprMean (now %g) or ' ...
        'MinCells (now %g).'], opts.MinExprMean, opts.MinCells);
end
X = X(keep, :);
T = T(keep, :);
if opts.Verbose
    fprintf('  %s passed the expression filter.\n', ...
        pkg.i_plural(height(T), 'gene'));
end

% ---- Library-size normalisation ----
libsize = full(sum(X, 1));
if any(libsize == 0)
    error(['Zero counts over the retained genes in %d of %d cells, which ' ...
        'cannot be normalised. Remove them first, e.g. with sc_qcfilter.'], ...
        sum(libsize == 0), nCells);
end
Y = full(X);
if opts.UseSingle
    Y = single(Y);
    libsize = single(libsize);
end
clear X
Y = (Y ./ libsize) * 10^round(log10(mean(double(libsize))));

% ---- Variance stabilisation, then log space ----
Y = 2 * sqrt(Y + 3/8);
Y = log2(Y + 1);

% ---- Clip at the mean per-cell extreme, symmetrically ----
threshold = mean(abs([mean(min(Y, [], 1)), mean(max(Y, [], 1))]));
Y = min(max(Y, -threshold), threshold);

% ---- Moving average along each chromosome ----
Y = i_smoothbychr(Y, T.Chr, opts.WindowLength, opts.Verbose);

% ---- Centre each cell, then express relative to the reference ----
Y = Y - median(Y, 1, 'omitnan');
refmeans = log2(mean(2.^Y(:, isref) - 1, 2) + 1);
Y = Y - refmeans;
Y = 2.^Y;

% ---- Flatten the range the reference cells themselves occupy ----
if opts.SdAmplifier > 0
    R = Y(:, isref);
    refmid = mean(R(:));
    refsd = mean(std(R, 0, 1)) * opts.SdAmplifier;
    clear R
    Y(Y > refmid - refsd & Y < refmid + refsd) = refmid;
end

% ---- Clamp residual outliers ----
Y = min(max(Y, mean(min(Y, [], 1))), mean(max(Y, [], 1)));

cnv = Y;

% ---- Optional piecewise-constant segmentation ----
seg = [];
if opts.Segment || nargout > 2
    if opts.Verbose
        fprintf('  segmenting...\n');
    end
    [segcnv, breakpoints] = pkg.e_cnvsegment(cnv, T.Chr, ...
        StopRatio=opts.SegStopRatio, UseSingle=opts.UseSingle);
    seg = struct('cnv', segcnv, 'breakpoints', breakpoints, ...
        'table', i_segtable(breakpoints, T));
end

if opts.Verbose
    fprintf('  done: cnv is %d genes by %d cells', size(cnv, 1), size(cnv, 2));
    if ~isempty(seg)
        fprintf(', in %s', pkg.i_plural(height(seg.table), 'segment'));
    end
    fprintf('.\n');
end
end


function S = i_segtable(breakpoints, T)
% One row per segment, in genome order. The last breakpoint is a sentinel
% one past the end, so segment k spans rows breakpoints(k) to
% breakpoints(k+1)-1.
nSeg = numel(breakpoints) - 1;
Chr = zeros(nSeg, 1);
FirstGene = strings(nSeg, 1);
LastGene = strings(nSeg, 1);
StartPos = zeros(nSeg, 1);
EndPos = zeros(nSeg, 1);
NumGenes = zeros(nSeg, 1);
for k = 1:nSeg
    rows = breakpoints(k):breakpoints(k+1) - 1;
    Chr(k) = T.Chr(rows(1));
    FirstGene(k) = T.Gene(rows(1));
    LastGene(k) = T.Gene(rows(end));
    StartPos(k) = T.Start(rows(1));
    EndPos(k) = T.Stop(rows(end));
    NumGenes(k) = numel(rows);
end
S = table(Chr, FirstGene, LastGene, StartPos, EndPos, NumGenes);
end


function isref = i_torefmask(refidx, nCells)
% Accept either a logical mask over cells or a list of column indices.
if islogical(refidx)
    if numel(refidx) ~= nCells
        error('REFIDX is a logical mask of length %d but X has %d columns.', ...
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


function Y = i_smoothbychr(Y, chr, windowLength, verbose)
% Triangular-kernel moving average, computed independently per chromosome so
% that no window ever straddles a centromere or a chromosome boundary. Edge
% genes are normalised by the kernel weight actually available to them
% rather than being padded with zeros.
if windowLength < 2
    warning('sc_infercnv:shortWindow', ...
        'WindowLength < 2 leaves the data unsmoothed.');
    return;
end
if mod(windowLength, 2) == 0
    windowLength = windowLength + 1;
    if verbose
        fprintf('  WindowLength must be odd; using %d.\n', windowLength);
    end
end

halfWindow = (windowLength - 1) / 2;
kernel = [1:halfWindow, halfWindow+1, halfWindow:-1:1]';
if isa(Y, 'single'), kernel = single(kernel); end

for c = unique(chr).'
    ix = find(chr == c);
    if numel(ix) < 2
        continue;
    end
    weight = conv2(ones(numel(ix), 1, 'like', Y), kernel, 'same');
    Y(ix, :) = conv2(Y(ix, :), kernel, 'same') ./ weight;
end
end
