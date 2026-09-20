function T_out = weight(T_ccc, X, genelist, c_cell_type, opts)
%GLY.WEIGHT  Re-weight cell-cell communication edges by glyco-state.
%   T_out = GLY.WEIGHT(T_ccc, X, genelist, c_cell_type) annotates an
%   existing ligand-receptor communication table with a glycosylation
%   modulation factor. Many protein LR interactions depend on the glycan
%   context of the participating cells (e.g. heparan sulfate co-receptors for
%   FGF/Wnt/chemokines, Fringe glycosylation of Notch, N-glycan branching for
%   EGFR-family residency). For each edge whose ligand/receptor family matches
%   the Channel-B map (GLY.LECTINMAP("lr")), the edge probability is
%   multiplied by a normalized per-cell-type score of the relevant glyco module
%   (GLY.STATE), taken from the compartment (sender or receiver) the
%   modulation acts on.
%
%   The operation is non-destructive: the original probability column is kept
%   and new columns are appended. Edges with no glyco dependency keep a factor
%   of 1.
%
%   T_out = GLY.WEIGHT(..., condition=COND, condLabels=["ctrl","treat"])
%   handles a two-condition table from SC_DOCK_CCC multi-group mode. Glyco
%   scores are then computed per cell type AND condition, each condition's
%   probability is weighted by its own factor, and the differential is
%   recomputed from the weighted values. Without this, a condition-averaged
%   factor would be applied to a between-condition difference, which cancels
%   out exactly the glycan remodelling the comparison is meant to detect.
%
%   USAGE:
%     % single-group
%     ccc = sc_dock_ccc(X, g, c);
%     T   = gly.weight(ccc.T_interactions, X, g, c);
%
%     % two-condition
%     ccc = sc_dock_ccc(X, g, c, "condition", cond, "cond_labels", {'ctrl','treat'});
%     T   = gly.weight(ccc.T_interactions, X, g, c, ...
%               condition=cond, condLabels=["ctrl","treat"]);
%
%   INPUTS:
%     T_ccc       - table with columns ligand, receptor, sender, receiver and a
%                   probability column (see probVar). One row per directed edge.
%     X           - genes-by-cells expression matrix (normalized)
%     genelist    - G-by-1 gene symbols (length = rows of X)
%     c_cell_type - N-by-1 cell-type label per cell (length = columns of X)
%
%   NAME-VALUE ARGUMENTS:
%     probVar    - name of the probability column in T_ccc (default "prob")
%     methodid   - scoring method for GLY.STATE (default 2)
%     minGenes   - minimum module genes present to score a module (default 3)
%     condition  - N-by-1 condition label per cell ([] = single-group mode)
%     condLabels - 1-by-2 condition names matching prob_cond1 / prob_cond2;
%                  required when condition is supplied and the table carries
%                  per-condition probabilities
%
%   OUTPUT:
%     T_out - copy of T_ccc with appended columns:
%             glyco_module  - modulating module ("" if none applied)
%             glyco_factor  - normalized [0,1] modulation factor (1 if none)
%             prob_weighted - <probVar> .* glyco_factor
%             ranked by prob_weighted descending.
%
%             In two-condition mode the appended columns are instead:
%             glyco_module, glyco_factor_cond1, glyco_factor_cond2,
%             prob_cond1_weighted, prob_cond2_weighted, delta_prob_weighted,
%             ranked by |delta_prob_weighted| descending.
%
%   The factor is min-max normalized per module across the groups present, so 0
%   marks the lowest-scoring group and 1 the highest. When several Channel-B
%   rows match an edge, the first match wins.
%
% see also: GLY.CCC, GLY.LECTINMAP, GLY.NORM,
%           GLY.LRWEIGHT, GLY.STATE, SC_DOCK_CCC

arguments
    T_ccc table
    X {mustBeNumeric}
    genelist (:, 1) string
    c_cell_type (:, 1) string
    opts.probVar (1, 1) string = "prob"
    opts.methodid (1, 1) double = 2
    opts.minGenes (1, 1) double = 3
    opts.condition (:, 1) string = strings(0, 1)
    opts.condLabels (:, 1) string = strings(0, 1)
end

vars = string(T_ccc.Properties.VariableNames);
required = ["ligand", "receptor", "sender", "receiver"];
missing = required(~ismember(required, vars));
if ~isempty(missing)
    error("GLY:WEIGHT:MissingColumns", ...
        "T_ccc is missing required column(s): %s.", strjoin(missing, ", "));
end
if numel(genelist) ~= size(X, 1) || numel(c_cell_type) ~= size(X, 2)
    error("GLY:WEIGHT:BadSizes", ...
        "GENELIST must match rows of X and C_CELL_TYPE must match columns of X.");
end

% Two-condition mode requires both a per-cell condition vector and a table that
% actually carries per-condition probabilities.
hasCondCols = all(ismember(["prob_cond1", "prob_cond2"], vars));
twoCond = ~isempty(opts.condition) && hasCondCols;

if ~isempty(opts.condition) && ~hasCondCols
    warning("GLY:WEIGHT:NoCondColumns", ...
        "CONDITION was supplied but T_ccc has no prob_cond1/prob_cond2 " + ...
        "columns; falling back to single-group weighting.");
end

% The mirror image of the guard above, which was missing. A table carrying
% prob_cond1/prob_cond2 came from a two-condition SC_DOCK_CCC run, and in
% that mode its probability column is abs(p2 - p1) -- a between-condition
% DIFFERENCE, not a pooled probability (sc_dock_ccc.m:423). Taking the
% single-group path then multiplies that difference by one
% condition-averaged glyco factor, which is the exact failure the header of
% this file warns about: a condition-averaged factor applied to a
% between-condition difference cancels out the glycan remodelling the
% comparison is meant to detect. It happened in silence, because TWOCOND
% needs both facts and only one of the two mismatches was reported.
% Measured on a two-row two-condition table with prob = 0.6 on both edges,
% omitting CONDITION returned prob_weighted = [0.6, 0] -- one interaction
% zeroed outright by an averaged factor.
if hasCondCols && isempty(opts.condition)
    warning("GLY:WEIGHT:CondColumnsUnused", ...
        "T_ccc has prob_cond1/prob_cond2 columns but CONDITION was not " + ...
        "supplied, so single-group weighting is used: one " + ...
        "condition-averaged factor is applied to the %s column, which " + ...
        "in a two-condition table is a between-condition difference. " + ...
        "Pass condition= and condLabels= for per-condition weighting.", ...
        opts.probVar);
end
if twoCond && numel(opts.condLabels) ~= 2
    error("GLY:WEIGHT:BadCondLabels", ...
        "CONDLABELS must name the two conditions matching prob_cond1 and " + ...
        "prob_cond2, e.g. condLabels=[""ctrl"",""treat""].");
end
if ~twoCond && ~ismember(opts.probVar, vars)
    error("GLY:WEIGHT:MissingColumns", ...
        "T_ccc is missing required column: %s.", opts.probVar);
end

% -------------------------------------------------------------------------
% Normalized glyco-module scores, grouped by cell type (and condition)
% -------------------------------------------------------------------------
if twoCond
    G = gly.norm(X, genelist, c_cell_type, condition=opts.condition, ...
        methodid=opts.methodid, minGenes=opts.minGenes);
else
    G = gly.norm(X, genelist, c_cell_type, ...
        methodid=opts.methodid, minGenes=opts.minGenes);
end

lig = string(T_ccc.ligand);
rec = string(T_ccc.receptor);
snd = string(T_ccc.sender);
rcv = string(T_ccc.receiver);

T_out = T_ccc;

% -------------------------------------------------------------------------
% Apply modulation per edge. minWeight=0 keeps the historical [0,1] factor
% range, in which a lowest-scoring compartment zeroes the edge out.
% -------------------------------------------------------------------------
if twoCond
    [w1, module] = gly.lrweight(lig, rec, snd, rcv, G, ...
        mode="raw", minWeight=0, condition=opts.condLabels(1));
    w2 = gly.lrweight(lig, rec, snd, rcv, G, ...
        mode="raw", minWeight=0, condition=opts.condLabels(2));

    T_out.glyco_module = module;
    T_out.glyco_factor_cond1 = w1;
    T_out.glyco_factor_cond2 = w2;
    T_out.prob_cond1_weighted = T_ccc.prob_cond1 .* w1;
    T_out.prob_cond2_weighted = T_ccc.prob_cond2 .* w2;
    T_out.delta_prob_weighted = T_out.prob_cond2_weighted - T_out.prob_cond1_weighted;

    [~, ord] = sort(abs(T_out.delta_prob_weighted), "descend");
    T_out = T_out(ord, :);
else
    [w, module] = gly.lrweight(lig, rec, snd, rcv, G, ...
        mode="raw", minWeight=0);

    T_out.glyco_module = module;
    T_out.glyco_factor = w;
    T_out.prob_weighted = T_ccc.(opts.probVar) .* w;
    T_out = sortrows(T_out, "prob_weighted", "descend");
end

end
