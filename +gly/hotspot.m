function [T, info] = hotspot(X, genelist, c_cell_type, opts)
%GLY.HOTSPOT  Which glycogenes are regulated between cell types, and where.
%
%   T = GLY.HOTSPOT(X, GENELIST, C_CELL_TYPE) splits the glycogene
%   repertoire into the part that every cell type runs at the same level and
%   the part that is regulated, and says for each regulated gene WHICH cell
%   types turn it up. Three metrics per gene, after Dworkin, Clausen & Joshi
%   (iScience 2022;25:104419, the Glycopacity paper):
%
%     ubiquity - the fraction of cell types expressing the gene at all.
%                Low ubiquity = the gene is switched on in a few cell types.
%     iqr      - the spread of its per-cell-type level. High spread = the
%                gene is present everywhere but run at different rates.
%     clr      - per (gene, cell type), how far above that gene's average
%                level across ALL cell types this one sits. This is the
%                metric that names the cell type, and the other two do not.
%
%   A gene is a hotspot if it is regulated by either route: restricted
%   expression (low ubiquity) or graded expression (high IQR / high CLR
%   somewhere). A gene that is neither is a "workhorse" - ubiquitously
%   expressed, stable, and therefore carrying no cell-type contrast at all.
%
%   USAGE:
%     T = gly.hotspot(sce.X, sce.g, sce.c_cell_type_tx);
%     [T, info] = gly.hotspot(X, g, ctype, MinCells=200);
%     T(T.class == "specific", :)          % on/off regulated glycogenes
%     info.clr(:, info.celltypes == "CAF") % what CAF turns up
%
%   WHY THE CLR MATTERS MORE THAN THE OTHER TWO HERE. Ubiquity and IQR are
%   per-GENE summaries across the whole dataset: they rank genes, and say
%   nothing about any one cell type. Every downstream consumer in this
%   toolbox that needs a glyco prior - GLY.NORM feeding GLY.LRWEIGHT and
%   TEN.SCTENIFOLDXCT_GLYCO - needs a per-(module, cell type) number, and
%   the CLR is the only one of the three that is one. It is also referenced
%   to the dataset's own mean rather than to the groups being compared,
%   which is what stops it collapsing when only two cell types are in play.
%
%   THE EXPRESSION CALL IS DETECTION-BASED, NOT THE PAPER'S CUT-OFF. The
%   paper calls a gene expressed by testing its pseudo-bulk value against a
%   cut-off of 0.0054, calibrated to 1 TPM on paired bulk/single-cell HEK293
%   data, with a cluster-size-dependent SD. That constant is only meaningful
%   at the end of their exact normalization chain, and recalibrating it is
%   not done here, so it is not used: a gene counts as expressed in a cell
%   type when it is detected in at least DetectFrac of that type's cells.
%   Detection fraction is in any case the better-behaved statistic for this
%   gene class (Chrysinas et al., NAR Genom Bioinform 2024;6:lqae169; see
%   GLY.ENRICH, which rests on the same argument). The consequence is that
%   `ubiquity` is comparable ACROSS CELL TYPES WITHIN a dataset but not
%   across datasets of different depth - check GLY.DETECT first, and do not
%   compare ubiquity between a 4600-UMI run and a 21000-UMI one.
%
%   CLR AND IQR ARE DEPTH-CORRECTED, UBIQUITY IS NOT. The per-cell-type
%   level is divided by that cell type's housekeeping level before either is
%   computed, which removes a per-cell-type offset that would otherwise add
%   itself to every gene's CLR in that column. Without it a deeply sequenced
%   cell type reads as turning up the entire glycogenome.
%
%   INPUTS:
%     X           - genes-by-cells RAW counts. The library normalization is
%                   done internally, so pre-normalized input double-counts it
%                   and makes the housekeeping correction a no-op.
%     genelist    - G-by-1 gene symbols (length = rows of X)
%     c_cell_type - N-by-1 cell-type label per cell (length = columns of X)
%
%   NAME-VALUE ARGUMENTS:
%     Collection - "genesets" (default, GLY.GENESETS) or "enzonto"
%                  (GLY.ENZONTO) for the glycogene universe and the module
%                  annotation attached to each gene.
%     Genes      - explicit glycogene list, overriding Collection's universe.
%     MinCells   - drop cell types with fewer cells (default 20). The paper
%                  uses 200, which is right for an atlas and often leaves
%                  nothing on a single experiment; raise it when you can.
%     DetectFrac - fraction of a cell type's cells in which a gene must be
%                  detected to count as expressed there (default 0.10).
%     ClrCut     - CLR above which a gene counts as high in a cell type
%                  (default 1, the paper's).
%     IqrCut     - IQR above which a ubiquitous gene counts as modulated
%                  (default 0.9, the paper's organ-level value).
%     SpecificCut- ubiquity at or below which a gene counts as cell-type
%                  specific (default 1/3, the paper's band edge; inclusive).
%     HKGenes    - housekeeping panel for the depth correction (default: see
%                  below). Only genes detected in a cell type contribute.
%     Trim       - percentile trim applied to the per-gene, per-cell-type
%                  mean (default [10 90], the paper's).
%
%   OUTPUTS:
%     T    - one row per glycogene, sorted by class then descending maxclr:
%            gene, module, nExpressed, ubiquity, iqr, maxclr, topCellType,
%            nHigh (cell types with CLR > ClrCut), class.
%     info - struct with the matrices behind T: .q (genes x cell types
%            pseudo-bulk), .qn (housekeeping-corrected), .clr, .detected
%            (fraction of cells), .expressed (logical), .genes, .celltypes,
%            .ncells, .hkfound, .hk (per-cell-type housekeeping level, the
%            depth offset that .qn removes), .opts.
%
%   CLASSES. "specific" (ubiquity <= SpecificCut), "modulated" (ubiquitous,
%   and IQR > IqrCut or CLR > ClrCut somewhere), "workhorse" (ubiquitous and
%   neither), "undetected" (expressed nowhere). The first two are hotspots.
%
%   THE SCALE IS COMPRESSED, SO ClrCut=1 IS A LARGE MOVE. The per-cell-type
%   level is a mean of log1p values after library normalization, so a raw
%   15-fold change in one gene's counts lands well under 1 on this scale.
%   That is the paper's scale and its cut-off, and on real data the range is
%   ample - GSE115978 gives maxclr up to 4.7 with 88 of 431 glycogenes
%   modulated - but a synthetic fixture with a modest planted contrast will
%   classify everything as a workhorse, which is the fixture being weak and
%   not the threshold being wrong.
%
%   THE DEFAULT HOUSEKEEPING PANEL IS NOT THE PAPER'S. Theirs is a 21-gene
%   panel shown only as a figure. This uses the Eisenberg & Levanon (Trends
%   Genet 2013;29:569-574) short list, which is the standard modern answer to
%   the same question, MINUS GPI - which is on their list and is also a
%   curated glycogene here (nucleotide sugar metabolism). Normalizing
%   glycogenes by a glycogene would partly cancel the signal being measured.
%
% see also: GLY.CAPACITY, GLY.DETECT, GLY.ENRICH, GLY.NORM, GLY.GENESETS

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    c_cell_type (:, 1) string
    opts.Collection (1, 1) string {mustBeMember(opts.Collection, ["genesets", "enzonto"])} = "genesets"
    opts.Genes (:, 1) string = strings(0, 1)
    opts.MinCells (1, 1) double {mustBePositive} = 20
    opts.DetectFrac (1, 1) double {mustBeInRange(opts.DetectFrac, 0, 1)} = 0.10
    opts.ClrCut (1, 1) double = 1
    opts.IqrCut (1, 1) double = 0.9
    opts.SpecificCut (1, 1) double {mustBeInRange(opts.SpecificCut, 0, 1)} = 1/3
    opts.HKGenes (:, 1) string = i_hkpanel()
    opts.Trim (1, 2) double = [10 90]
end

if numel(genelist) ~= size(X, 1)
    error("GLY:HOTSPOT:BadGenelist", ...
        "GENELIST length (%d) must equal the number of rows of X (%d).", ...
        numel(genelist), size(X, 1));
end
if numel(c_cell_type) ~= size(X, 2)
    error("GLY:HOTSPOT:BadLabels", ...
        "C_CELL_TYPE length (%d) must equal the number of columns of X (%d).", ...
        numel(c_cell_type), size(X, 2));
end
if opts.Trim(1) >= opts.Trim(2) || opts.Trim(1) < 0 || opts.Trim(2) > 100
    error("GLY:HOTSPOT:BadTrim", "TRIM must be [lo hi] with 0 <= lo < hi <= 100.");
end

% -------------------------------------------------------------------------
% Cell types that clear MinCells
% -------------------------------------------------------------------------
[types, ~, typeIdx] = unique(c_cell_type);
ncells = accumarray(typeIdx, 1);
keepType = ncells >= opts.MinCells;
if ~any(keepType)
    error("GLY:HOTSPOT:NoCellTypes", ...
        "No cell type has at least MinCells=%d cells (largest has %d).", ...
        opts.MinCells, max(ncells));
end
if sum(keepType) < 2
    error("GLY:HOTSPOT:OneCellType", ...
        ['Only one cell type clears MinCells=%d. Every metric here is a ' ...
        'contrast BETWEEN cell types and none is defined on one.'], opts.MinCells);
end
types = types(keepType);
ncells = ncells(keepType);

% -------------------------------------------------------------------------
% The glycogene universe, and each gene's module annotation
% -------------------------------------------------------------------------
[glyGenes, geneModule] = i_universe(opts);

% Match into the data, case-insensitively (mouse orthologues share spelling)
upGL = upper(genelist);
[inData, rowOf] = ismember(upper(glyGenes), upGL);
if ~any(inData)
    error("GLY:HOTSPOT:NoGlycogenes", ...
        "None of the %d glycogenes in the %s collection are in GENELIST.", ...
        numel(glyGenes), opts.Collection);
end
glyGenes = glyGenes(inData);
geneModule = geneModule(inData);
glyRows = rowOf(inData);

% Housekeeping rows (panel genes present in the data)
[hkIn, hkRow] = ismember(upper(opts.HKGenes), upGL);
hkfound = opts.HKGenes(hkIn);
hkRows = hkRow(hkIn);
if numel(hkfound) < 3
    error("GLY:HOTSPOT:NoHousekeeping", ...
        ['Only %d of the %d housekeeping genes are in GENELIST. The depth ' ...
        'correction is not trustworthy below 3; pass HKGenes explicitly.'], ...
        numel(hkfound), numel(opts.HKGenes));
end

% -------------------------------------------------------------------------
% Per-cell library normalization to 10k counts, then log1p
% -------------------------------------------------------------------------
% Only the glycogene and housekeeping rows are ever read, but the library
% size must come from the FULL matrix - a depth computed over glycogenes
% alone would be the very quantity being corrected for.
libsize = full(sum(X, 1));
libsize(libsize == 0) = 1;
scale = 1e4 ./ libsize;

Xg = full(X(glyRows, :)) .* scale;
Xg = log1p(Xg);
Xh = full(X(hkRows, :)) .* scale;
Xh = log1p(Xh);

% -------------------------------------------------------------------------
% Pseudo-bulk per cell type: trimmed mean over DETECTED cells only
% -------------------------------------------------------------------------
nGly = numel(glyGenes);
nType = numel(types);
q = zeros(nGly, nType);
detected = zeros(nGly, nType);
qhk = zeros(numel(hkfound), nType);

typeOf = c_cell_type;
for k = 1:nType
    sel = typeOf == types(k);
    q(:, k) = i_trimmean(Xg(:, sel), opts.Trim);
    detected(:, k) = mean(Xg(:, sel) > 0, 2);
    qhk(:, k) = i_trimmean(Xh(:, sel), opts.Trim);
end

% -------------------------------------------------------------------------
% Housekeeping correction: a per-cell-type offset, removed on the log scale
% -------------------------------------------------------------------------
% Geometric mean over the panel genes DETECTED in that cell type, as in the
% paper - the point is robustness to a panel gene dropping out of one column,
% which would otherwise move every glycogene in that column.
h = zeros(1, nType);
for k = 1:nType
    v = qhk(:, k);
    v = v(v > 0);
    if isempty(v)
        error("GLY:HOTSPOT:DeadHousekeeping", ...
            "No housekeeping gene is detected in cell type '%s'.", types(k));
    end
    h(k) = exp(mean(log(v)));
end

qn = nan(nGly, nType);
for k = 1:nType
    idx = q(:, k) > 0;
    qn(idx, k) = log2(q(idx, k) ./ h(k));
end

% -------------------------------------------------------------------------
% The three metrics
% -------------------------------------------------------------------------
expressed = detected >= opts.DetectFrac;
ubiquity = mean(expressed, 2);

% CLR: centre each gene's log level across the cell types where it is
% measurable. A gene absent from a column has no level there, and imputing
% one (zero, or a floor) would invent a depression the data did not show.
clr = nan(nGly, nType);
iqrv = nan(nGly, 1);
for i = 1:nGly
    v = qn(i, :);
    ok = ~isnan(v);
    if nnz(ok) < 2
        continue;
    end
    clr(i, ok) = v(ok) - mean(v(ok));
    iqrv(i) = iqr(v(ok));
end

nHigh = sum(clr > opts.ClrCut, 2);
[maxclr, topIdx] = max(clr, [], 2, "omitnan");
topCellType = types(topIdx);
allNaN = all(isnan(clr), 2);
maxclr(allNaN) = NaN;

% A gene measurable in only ONE cell type has no CLR - centring needs two
% columns - but it is the most cell-type-specific case there is, and
% returning no cell type for it would blank exactly the rows a caller came
% for. Name the cell type from the detection instead, and leave maxclr NaN
% so nobody mistakes it for a measured contrast.
if any(allNaN)
    [bestDet, detIdx] = max(detected(allNaN, :), [], 2);
    fallback = types(detIdx);
    fallback(bestDet == 0) = "";        % detected nowhere: nothing to name
    topCellType(allNaN) = fallback;
end

% -------------------------------------------------------------------------
% Classification
% -------------------------------------------------------------------------
class = repmat("workhorse", nGly, 1);
class(ubiquity == 0) = "undetected";
% Inclusive at the band edge. Ubiquity can only take the values j/K, so a
% strict < makes the boundary behave erratically as K changes: at K=3 no
% gene can EVER be specific (1/3 is not < 1/3), at K=6 it takes one cell
% type, at K=9 three. Inclusive gives the monotone reading - at most a third
% of the cell types - at every K.
isSpecific = ubiquity > 0 & ubiquity <= opts.SpecificCut;
class(isSpecific) = "specific";
isModulated = ~isSpecific & ubiquity > 0 & ...
    ((iqrv > opts.IqrCut) | (nHigh > 0));
class(isModulated) = "modulated";

T = table(glyGenes(:), geneModule(:), sum(expressed, 2), ubiquity, ...
    iqrv, maxclr, topCellType(:), nHigh, class, ...
    VariableNames = ["gene", "module", "nExpressed", "ubiquity", ...
    "iqr", "maxclr", "topCellType", "nHigh", "class"]);

classOrder = ["specific", "modulated", "workhorse", "undetected"];
[~, ord] = ismember(T.class, classOrder);
T.sortkey = ord;
T = sortrows(T, ["sortkey", "maxclr"], ["ascend", "descend"]);
T.sortkey = [];

info = struct();
info.q = q;
info.qn = qn;
info.clr = clr;
info.detected = detected;
info.expressed = expressed;
info.genes = glyGenes(:);
info.modules = geneModule(:);
info.celltypes = types(:);
info.ncells = ncells(:);
info.hkfound = hkfound(:);
info.hk = h(:);          % per-cell-type housekeeping level, the depth offset
info.opts = opts;

end


% =========================================================================
function m = i_trimmean(V, trim)
% Per-row mean over detected (nonzero) values, dropping those outside the
% trim percentiles of that row's detected values. Rows with no detected
% value return 0.
n = size(V, 1);
m = zeros(n, 1);
for i = 1:n
    v = V(i, :);
    v = v(v > 0);
    if isempty(v)
        continue;
    end
    if numel(v) > 2
        lo = prctile(v, trim(1));
        hi = prctile(v, trim(2));
        vt = v(v >= lo & v <= hi);
        if ~isempty(vt)
            v = vt;
        end
    end
    m(i) = mean(v);
end
end


% =========================================================================
function [genes, modules] = i_universe(opts)
% The glycogene list, plus one module label per gene. A gene in several
% modules is labelled with all of them, joined by ";".
if ~isempty(opts.Genes)
    genes = unique(upper(opts.Genes));
    modules = repmat("", numel(genes), 1);
    return;
end

if opts.Collection == "enzonto"
    [setmatrx, setnames, setgenes] = gly.enzonto();
else
    [setmatrx, setnames, setgenes] = gly.genesets();
end

genes = string(setgenes(:));
modules = strings(numel(genes), 1);
setnames = string(setnames(:));
for i = 1:numel(genes)
    mem = setnames(logical(setmatrx(:, i)));
    if ~isempty(mem)
        modules(i) = join(mem, ";");
    end
end
end


% =========================================================================
function g = i_hkpanel()
% Eisenberg & Levanon (Trends Genet 2013;29:569-574) short list, less GPI,
% which is itself a curated glycogene here. See the header note.
g = ["C1orf43"; "CHMP2A"; "EMC7"; "PSMB2"; "PSMB4"; "RAB7A"; ...
    "REEP5"; "SNRPD3"; "VCP"; "VPS29"];
end
