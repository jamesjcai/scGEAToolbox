function T_out = shield(receptorGenes, X, genelist, c_cell_type, opts)
%GLY.SHIELD  Cell-type-specific N-glycan shielding of docking targets.
%   T = GLY.SHIELD(RECEPTORGENES, X, GENELIST, C_CELL_TYPE) estimates how
%   strongly each receptor's surface is shielded by N-linked glycans in each
%   cell type, and optionally discounts molecular-docking affinities by that
%   amount.
%
%   Docking is run against PDB structures, which are almost always solved
%   deglycosylated or have their glycans stripped during preparation. A
%   compound can therefore score well against a bare structure while never
%   reaching the receptor on a real cell whose glycocalyx covers the site. The
%   shielding score combines two independent pieces of evidence:
%
%     1. How glycosylatable the protein is - the count of N-X-S/T sequons
%        (X not proline) in its UniProt sequence, the canonical N-glycosylation
%        motif.
%     2. How much N-glycosylation machinery the cell type actually runs - the
%        normalized N-glycan module scores from GLY.NORM.
%
%   A protein full of sequons in a cell type with no N-glycan biosynthesis is
%   not shielded, and neither is a sequon-free protein in a heavily
%   glycosylating cell. Only the product is informative.
%
%   USAGE:
%     run.ml_scDock();
%     ccc = sc_dock_ccc(X, g, ctype);
%     T   = gly.shield(unique(ccc.T_interactions.receptor), X, g, ctype);
%
%     % discount docking affinities for one cell type
%     T = gly.shield(genes, X, g, ctype, cellType="Tumor", affinity=aff);
%
%   INPUTS:
%     receptorGenes - R-by-1 receptor gene symbols to evaluate
%     X             - genes-by-cells expression matrix (normalized)
%     genelist      - G-by-1 gene symbols (length = rows of X)
%     c_cell_type   - N-by-1 cell-type label per cell (length = columns of X)
%
%   NAME-VALUE ARGUMENTS:
%     cellType  - restrict to these cell types (default: all)
%     alpha     - strength of the affinity discount, in [0,1] (default 0.5)
%     affinity  - R-by-1 docking affinities in kcal/mol aligned with
%                 RECEPTORGENES; when given, affinity_adj is appended
%     sequences - R-by-1 protein sequences aligned with RECEPTORGENES, to skip
%                 the UniProt lookup (useful offline or for non-human data)
%     modules   - N-glycan module names to average for cellular capacity
%                 (default: the ER biosynthesis and Golgi processing modules)
%     saturation- sequon count at which shielding reaches 1-1/e (default 5)
%
%   OUTPUT:
%     T_out - table with one row per receptor x cell type:
%       receptor, n_sequons, seq_length, sequon_density, shield_capacity,
%       cell_type, glyco_capacity, shield_score [, affinity, affinity_adj]
%     ranked by shield_score descending.
%
%   Shielding from sequon count uses a saturating absolute transform,
%   1 - exp(-n_sequons/saturation), rather than normalizing across the queried
%   genes. The score for a given protein therefore does not change when other
%   genes are added to or removed from the query, which matters if the values
%   are to be reported or compared across analyses.
%
%   The adjusted affinity is affinity*(1 - alpha*shield_score). Vina affinities
%   are negative and more negative means tighter, so shielding moves the value
%   toward zero, weakening the predicted interaction.
%
%   This is a surface-level rather than pocket-level estimate: SC_DOCK_VINA
%   docks blind over a whole-protein bounding box, so there is no defined
%   pocket to measure sequon distance from. Global shielding is the matched
%   granularity for that search. If a pocket is ever specified, weighting
%   sequons by distance to it would be the natural refinement.
%
%   Sequences are fetched once per gene per session and cached.
%
% see also: SC_DOCK_VINA, SC_DOCK_GENE2PDB, SC_DOCK_CCC, GLY.NORM,
%           GLY.STATE, RUN.ML_SCDOCK

arguments
    receptorGenes (:, 1) string
    X {mustBeNumeric}
    genelist (:, 1) string
    c_cell_type (:, 1) string
    opts.cellType (:, 1) string = strings(0, 1)
    opts.alpha (1, 1) double {mustBeBetween(opts.alpha, 0, 1)} = 0.5
    opts.affinity (:, 1) double = []
    opts.sequences (:, 1) string = strings(0, 1)
    opts.modules (:, 1) string = ["Glyco_N_glycan_biosynthesis_ER"; ...
                                  "Glyco_N_glycan_processing_Golgi"]
    opts.saturation (1, 1) double {mustBePositive} = 5
end

nGene = numel(receptorGenes);
if ~isempty(opts.affinity) && numel(opts.affinity) ~= nGene
    error("GLY:SHIELD:BadAffinity", ...
        "AFFINITY must have one value per receptor gene (%d), has %d.", ...
        nGene, numel(opts.affinity));
end
if ~isempty(opts.sequences) && numel(opts.sequences) ~= nGene
    error("GLY:SHIELD:BadSequences", ...
        "SEQUENCES must have one entry per receptor gene (%d), has %d.", ...
        nGene, numel(opts.sequences));
end

% -------------------------------------------------------------------------
% Sequence-side: N-X-S/T sequon counts
% -------------------------------------------------------------------------
nSequon = nan(nGene, 1);
seqLen = nan(nGene, 1);
for k = 1:nGene
    if ~isempty(opts.sequences)
        seq = opts.sequences(k);
    else
        seq = i_fetchsequence(receptorGenes(k));
    end
    if strlength(seq) == 0
        continue;       % leave as NaN; reported but not scored
    end
    seqLen(k) = strlength(seq);
    nSequon(k) = i_countsequons(seq);
end

missing = isnan(seqLen);
if any(missing)
    warning("GLY:SHIELD:NoSequence", ...
        "No sequence retrieved for %d gene(s): %s", sum(missing), ...
        strjoin(receptorGenes(missing), ", "));
end

sequonDensity = nSequon ./ seqLen;
shieldCapacity = 1 - exp(-nSequon./opts.saturation);

% -------------------------------------------------------------------------
% Cell-side: N-glycan biosynthesis capacity per cell type
% -------------------------------------------------------------------------
G = gly.norm(X, genelist, c_cell_type);

modRows = find(ismember(G.setnames, opts.modules));
if isempty(modRows)
    error("GLY:SHIELD:NoModules", ...
        "None of the requested N-glycan modules were scored: %s.", ...
        strjoin(opts.modules, ", "));
end
capacity = mean(G.scores(modRows, :), 1);      % 1-by-nGroups

cellTypes = G.celltypes;
if ~isempty(opts.cellType)
    keep = ismember(cellTypes, opts.cellType);
    if ~any(keep)
        error("GLY:SHIELD:NoSuchCellType", ...
            "None of the requested cell types are present: %s.", ...
            strjoin(opts.cellType, ", "));
    end
    cellTypes = cellTypes(keep);
    capacity = capacity(keep);
end
nCt = numel(cellTypes);

% -------------------------------------------------------------------------
% Cross receptors with cell types
% -------------------------------------------------------------------------
[gi, ci] = ndgrid(1:nGene, 1:nCt);
gi = gi(:);
ci = ci(:);

T_out = table(receptorGenes(gi), nSequon(gi), seqLen(gi), sequonDensity(gi), ...
    shieldCapacity(gi), cellTypes(ci), capacity(ci)', ...
    shieldCapacity(gi).*capacity(ci)', ...
    VariableNames=["receptor", "n_sequons", "seq_length", "sequon_density", ...
    "shield_capacity", "cell_type", "glyco_capacity", "shield_score"]);

if ~isempty(opts.affinity)
    aff = opts.affinity(gi);
    T_out.affinity = aff;
    T_out.affinity_adj = aff.*(1 - opts.alpha.*T_out.shield_score);
end

T_out = sortrows(T_out, "shield_score", "descend");

end


%% ---- count N-X-S/T sequons (X not proline), including overlaps ----
function n = i_countsequons(seq)
% Lookahead keeps overlapping motifs such as NNSS from being missed.
hits = regexp(upper(char(seq)), 'N(?=[^P][ST])');
n = numel(hits);
end


%% ---- UniProt sequence lookup, cached per session ----
function seq = i_fetchsequence(gene)
persistent cache
if isempty(cache)
    cache = containers.Map("KeyType", "char", "ValueType", "char");
end
key = char(upper(gene));
if isKey(cache, key)
    seq = string(cache(key));
    return;
end

seq = "";
try
    url = sprintf(['https://rest.uniprot.org/uniprotkb/search?' ...
        'query=gene:%s+organism_id:9606+reviewed:true' ...
        '&format=tsv&fields=accession,sequence'], urlencode(key));
    txt = webread(url, weboptions("Timeout", 30));
    lines = strsplit(strtrim(txt), newline);
    if numel(lines) >= 2
        cols = strsplit(lines{2}, char(9));
        if numel(cols) >= 2
            seq = string(strtrim(cols{2}));
        end
    end
catch
    % Leave seq empty; the caller reports the gene as unresolved.
end

cache(key) = char(seq);

end
