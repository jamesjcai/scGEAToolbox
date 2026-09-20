function G = norm(X, genelist, c_cell_type, opts)
%NORM  Normalized per-group glycobiological state scores.
%   G = GLY.NORM(X, GENELIST, C_CELL_TYPE) scores every cell against the
%   curated glycobiology modules (GLY.STATE), averages the scores within
%   each cell type, and min-max normalizes each module across the groups. The
%   result is a compact, reusable readout of "how much of module M can group K
%   make, relative to the other groups in this dataset".
%
%   G = GLY.NORM(..., condition=COND) groups by cell type AND condition
%   instead of cell type alone. Use this whenever the downstream analysis
%   compares two conditions: glycosylation programs are themselves remodelled
%   between conditions, so a condition-averaged factor would cancel out the
%   very signal being tested.
%
%   USAGE:
%     G = gly.norm(X, sce.g, sce.c_cell_type_tx);
%     G = gly.norm(X, sce.g, sce.c_cell_type_tx, condition=cond);
%
%   INPUTS:
%     X           - genes-by-cells expression matrix (normalized, e.g. log1p)
%     genelist    - G-by-1 gene symbols (length = rows of X)
%     c_cell_type - N-by-1 cell-type label per cell (length = columns of X)
%
%   NAME-VALUE ARGUMENTS:
%     condition - N-by-1 condition label per cell; [] = group by cell type only
%     methodid  - scoring method for GLY.STATE (1 UCell, 2 AddModuleScore
%                 default, 3 AUCell)
%     minGenes  - minimum module genes present to score a module (default 3)
%
%   OUTPUT:
%     G - struct with fields:
%       .scores     - M-by-K normalized scores in [0,1] (M modules x K groups)
%       .raw        - M-by-K unnormalized group means
%       .setnames   - M-by-1 module names
%       .groups     - K-by-1 group keys ("celltype" or "celltype|condition")
%       .celltypes  - K-by-1 cell type of each group
%       .conditions - K-by-1 condition of each group ("" when not used)
%       .hascond    - logical, true when grouping included condition
%
%   Normalization is per module across groups, so 0 marks the lowest-scoring
%   group and 1 the highest. A module with no spread across groups maps to 0
%   for every group; callers that need a neutral default should treat an
%   all-zero module as "uninformative" rather than "absent".
%
%   Use GLY.GROUPCOL to look up the column for a given cell type and
%   condition.
%
%   NORMALIZE="CLR": KEEP THE EFFECT SIZE INSTEAD OF DISCARDING IT.
%
%   With min-max, the lowest group is 0 and the highest is 1 BY
%   CONSTRUCTION, whatever the gap between them. On two groups that is the
%   whole output: every module reads exactly 0 and 1, and a module whose two
%   group means differ by 0.001 is indistinguishable from one where they
%   differ by 10. The magnitude is not attenuated, it is erased, and every
%   downstream weight is then a function of the sign of the difference only.
%
%   NORMALIZE="clr" replaces this with the centred log ratio of Dworkin,
%   Clausen & Joshi (iScience 2022;25:104419): each module's per-group level
%   is divided by that group's housekeeping level (removing a per-group depth
%   offset), logged, and centred across the groups. A logistic squash then
%   puts it back on [0,1] for GLY.LRWEIGHT, with 0.5 - and hence a
%   modulation factor of exactly 1 - at the module's dataset-wide mean.
%
%   Three consequences worth stating plainly:
%
%     - 0.5 means "average for this module in this dataset", a defensible
%       neutral. Min-max's 0 means "lowest of the groups present", which is
%       not a neutral at all and drifts with which groups were passed in.
%     - The scale is preserved: two groups give symmetric scores around 0.5
%       whose distance from it tracks the actual log fold difference.
%     - The reference is every group in the call. Pass the WHOLE dataset and
%       read two columns off it, rather than subsetting to two cell types
%       first - subsetting first makes the reference a mean of two and
%       throws away the advantage, though the effect-size one survives.
%
%   The two modes also differ in the underlying quantity, not just its
%   rescaling: "minmax" averages per-cell AddModuleScores (GLY.STATE),
%   "clr" aggregates a per-group pseudo-bulk (GLY.HOTSPOT). AddModuleScore
%   subtracts a dataset-matched control set and is pinned near zero by
%   construction, which is why it is not a sound input to a log ratio.
%
%   In "clr" mode G.raw holds the centred log2 ratios themselves, so
%   G.raw = 0 marks the module's dataset mean and +1 a doubling above it.
%
% see also: GLY.STATE, GLY.HOTSPOT, GLY.GROUPCOL, GLY.LRWEIGHT,
%           GLY.WEIGHT, TEN.SCTENIFOLDXCT

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    c_cell_type (:, 1) string
    opts.condition (:, 1) string = strings(0, 1)
    opts.methodid (1, 1) double = 2
    opts.minGenes (1, 1) double = 3
    opts.normalize (1, 1) string {mustBeMember(opts.normalize, ["minmax", "clr"])} = "minmax"
    opts.clrScale (1, 1) double {mustBePositive} = 1
end

if numel(genelist) ~= size(X, 1)
    error("GLY:NORM:BadGenelist", ...
        "GENELIST length (%d) must equal the number of rows of X (%d).", ...
        numel(genelist), size(X, 1));
end
if numel(c_cell_type) ~= size(X, 2)
    error("GLY:NORM:BadLabels", ...
        "C_CELL_TYPE length (%d) must equal the number of columns of X (%d).", ...
        numel(c_cell_type), size(X, 2));
end

hascond = ~isempty(opts.condition);
if hascond && numel(opts.condition) ~= size(X, 2)
    error("GLY:NORM:BadCondition", ...
        "CONDITION length (%d) must equal the number of columns of X (%d).", ...
        numel(opts.condition), size(X, 2));
end

if hascond
    groupkey = c_cell_type + "|" + opts.condition;
else
    groupkey = c_cell_type;
end

if opts.normalize == "clr"
    [scores, raw, setnames, groups] = i_clrscores(X, genelist, groupkey, opts);
else
    % ---------------------------------------------------------------------
    % Per-cell module scores, then group means
    % ---------------------------------------------------------------------
    [cs, setnames] = gly.state(X, genelist, opts.methodid, opts.minGenes);

    raw = pkg.i_grpmean(cs, groupkey);   % modules x groups (groups sorted)
    groups = unique(groupkey);           % matches i_grpmean column order

    % ---------------------------------------------------------------------
    % Min-max normalize each module across groups
    % ---------------------------------------------------------------------
    lo = min(raw, [], 2);
    hi = max(raw, [], 2);
    span = hi - lo;
    span(span == 0) = 1;                 % avoid divide-by-zero
    scores = (raw - lo) ./ span;
end

% -------------------------------------------------------------------------
% Split the group keys back into their parts
% -------------------------------------------------------------------------
if hascond
    parts = split(groups, "|");
    celltypes = parts(:, 1);
    conditions = parts(:, 2);
else
    celltypes = groups;
    conditions = repmat("", numel(groups), 1);
end

G.scores = scores;
G.raw = raw;
G.setnames = setnames(:);
G.groups = groups(:);
G.celltypes = celltypes(:);
G.conditions = conditions(:);
G.hascond = hascond;
G.normalize = opts.normalize;

end


% =========================================================================
function [scores, clr, setnames, groups] = i_clrscores(X, genelist, groupkey, opts)
%I_CLRSCORES  Module-level centred log ratio across groups, squashed to [0,1].
%
% The per-group, per-gene pseudo-bulk and the per-group housekeeping level
% both come from GLY.HOTSPOT, so the two functions cannot drift apart in how
% they define a group's expression level.

[setmatrx, allnames, setgenes] = gly.genesets();
setnames = string(allnames(:));

% MinCells=1: GLY.NORM has never dropped a group, and callers index its
% output by cell type through GLY.GROUPCOL. Dropping groups here would make
% that lookup silently return 0 for a small cell type.
[~, H] = gly.hotspot(X, genelist, groupkey, ...
    Genes = string(setgenes(:)), MinCells = 1);

groups = H.celltypes(:);
nGroup = numel(groups);

% Map each module onto the rows of H.q that its genes occupy
hGenes = upper(string(H.genes(:)));
uGenes = upper(string(setgenes(:)));

nSet = numel(setnames);
level = nan(nSet, nGroup);
keep = false(nSet, 1);
for m = 1:nSet
    gm = uGenes(logical(setmatrx(m, :)));
    [inH, row] = ismember(gm, hGenes);
    rows = row(inH);
    if numel(rows) < opts.minGenes
        continue;       % same guard as GLY.STATE: too few genes to score
    end
    keep(m) = true;
    % Level WHERE EXPRESSED times FRACTION EXPRESSING, i.e. the
    % unconditional group mean, not the pseudo-bulk on its own.
    %
    % This is not a detail. GLY.HOTSPOT's q is a trimmed mean over the cells
    % that detect the gene, so it says how high a gene runs GIVEN that it
    % runs at all, and averaging that across a module gives a cell type full
    % credit for a gene seen in 2% of its cells. Measured on GSE115978, the
    % conditional version put T.CD4 ABOVE CAF on the heparan-sulfate module
    % - contradicting GLY.CAPACITY, which finds HS closed in T.CD4 (no EXT1)
    % and open in CAF, on the same data. Multiplying by the detection
    % fraction restores CAF 1.64 > Endo. 1.33 > Mal 0.71 > every lymphocyte,
    % which is both the known biology and the capacity call. Detection
    % fraction is in any case the better-behaved statistic for this gene
    % class - the argument GLY.ENRICH rests on.
    level(m, :) = mean(H.q(rows, :) .* H.detected(rows, :), 1);
end

if ~any(keep)
    error("GLY:NORM:NoModules", ...
        "No glyco module has at least minGenes=%d genes present in the data.", ...
        opts.minGenes);
end
setnames = setnames(keep);
level = level(keep, :);

% Housekeeping correction, then log. The floor keeps a module that is
% entirely undetected in one group finite; it sits well below any real
% level (the pseudo-bulk of a detected gene is a mean of log1p values).
floorq = 1e-6;
hk = H.hk(:)';
logLevel = log2((level + floorq) ./ hk);

% Centre across groups: this is the CLR proper.
clr = logLevel - mean(logLevel, 2);

% Logistic squash back onto [0,1] for GLY.LRWEIGHT, 0.5 at the dataset mean.
scores = 1 ./ (1 + exp(-clr ./ opts.clrScale));

end
