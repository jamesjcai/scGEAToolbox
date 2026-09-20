function [w, moduleName, compartment, detail] = lrweight(ligNames, recNames, senderTypes, receiverTypes, G, opts)
%LRWEIGHT  Glyco-modulation weight for protein ligand-receptor pairs.
%   [W, MODULE] = GLY.LRWEIGHT(LIG, REC, SENDER, RECEIVER, G) returns a
%   per-pair multiplicative weight reflecting whether the glycan context needed
%   by each ligand-receptor interaction is actually available in the relevant
%   compartment. Pairs matching a Channel-B row of GLY.LECTINMAP (e.g.
%   FGF-FGFR requiring heparan sulfate on the receiver) are scaled by that
%   compartment's normalized module score; all other pairs get a weight of 1.
%
%   USAGE:
%     G = gly.norm(X, g, ctype);
%     w = gly.lrweight(T.ligand, T.receptor, T.sender, T.receiver, G);
%
%   INPUTS:
%     ligNames      - N-by-1 ligand gene symbols
%     recNames      - N-by-1 receptor gene symbols
%     senderTypes   - N-by-1 sender cell types, or a scalar to broadcast
%     receiverTypes - N-by-1 receiver cell types, or a scalar to broadcast
%     G             - struct from GLY.NORM
%
%   NAME-VALUE ARGUMENTS:
%     lambda    - modulation strength (default 0.5)
%     mode      - "center" (default) or "boost"; see below
%     condition - scalar condition label used to select the group column when
%                 G was built with conditions (default "")
%     minWeight - floor applied to the returned weight (default 0.05)
%     usenode   - apply the node-level factor below (default true)
%
%   OUTPUTS:
%     w           - N-by-1 weights
%     moduleName  - N-by-1 modulating module name ("" where no rule applied)
%     compartment - N-by-1 "sender"/"receiver" ("" where no rule applied)
%
%   Three weighting modes, differing in what counts as neutral:
%     "center" : w = 1 + lambda*(2f - 1), so the lowest-scoring cell type is
%                down-weighted (1-lambda) and the highest is up-weighted
%                (1+lambda). Use this when the question is whether a glycan
%                requirement is met, since an unmet requirement should cost the
%                pair something relative to pairs with no requirement at all.
%     "boost"  : w = 1 + lambda*f, so no pair is ever penalized below the
%                unmodulated baseline of 1. Use this when the glyco annotation
%                is treated as supporting evidence only.
%     "raw"    : w = f, ignoring lambda. This is the [0,1] attenuation factor
%                used by GLY.WEIGHT, where a matched pair is scaled down in
%                proportion to the compartment's capacity.
%
%   Because module scores from GLY.NORM are min-max normalized across
%   the groups present, f is a relative rather than absolute capacity: f=0 means
%   "lowest among the cell types in this dataset", not "no biosynthetic
%   capacity". Interpret "center" mode accordingly.
%
%   NODE-LEVEL FACTOR (phi_node), MULTIPLIED INTO THE SAME phi. Rules match by
%   symbol PATTERN (e.g. "^FGFR"), so every family member a rule matches has
%   always received the identical module-level factor above - FGFR1 and FGFR4
%   were indistinguishable to Channel B, and on GSE115978 this was 209 of 255
%   tested pairs sharing one value (DOCS/RESULTS_GSE115978.md). phi_node breaks
%   that tie using each gene's OWN documented glycosylation site count
%   (GLY.NODEANNOT, from UniProt): the gene read is the ligand when the
%   rule's compartment is "sender", the receptor when it is "receiver" - the
%   same convention MAPB.COMPARTMENT already uses to pick which cell population
%   the module score comes from, reused here rather than newly curated per
%   rule. Site counts are min-max normalized ACROSS THE SET OF DISTINCT GENES
%   THE CURRENT RULE MATCHED IN THIS CALL (not globally), so phi_node=0.5 (->
%   factor 1 in "center" mode) whenever a rule matches only one distinct gene
%   or every matched gene happens to carry the same count - there is no
%   information to discriminate on and none is invented. The same
%   MODE/LAMBDA/EXPONENT apply, so lambda=0 collapses phi_node to 1 exactly
%   like phi_module, and the lambda=0 control identity is unaffected. Set
%   usenode=false to disable and reproduce the pre-node-annotation behaviour.
%
% see also: GLY.NORM, GLY.LECTINMAP, GLY.NODEANNOT,
%           GLY.WEIGHT, TEN.SCTENIFOLDXCT

arguments
    ligNames (:, 1) string
    recNames (:, 1) string
    senderTypes (:, 1) string
    receiverTypes (:, 1) string
    G (1, 1) struct
    opts.lambda (1, 1) double = 0.5
    opts.mode (1, 1) string {mustBeMember(opts.mode, ["center", "boost", "raw"])} = "center"
    opts.condition (1, 1) string = ""
    opts.minWeight (1, 1) double = 0.05
    opts.usenode (1, 1) logical = true
end

n = numel(ligNames);
if numel(recNames) ~= n
    error("GLY:LRWEIGHT:BadSizes", ...
        "LIGNAMES (%d) and RECNAMES (%d) must have the same length.", ...
        n, numel(recNames));
end
senderTypes = i_broadcast(senderTypes, n, "SENDERTYPES");
receiverTypes = i_broadcast(receiverTypes, n, "RECEIVERTYPES");

w = ones(n, 1);
moduleName = repmat("", n, 1);
compartment = repmat("", n, 1);
if n == 0
    return;
end

mapB = gly.lectinmap("lr");

if opts.usenode
    nodeT = gly.nodeannot();
end

% -------------------------------------------------------------------------
% Resolve each distinct cell type to a group column once
% -------------------------------------------------------------------------
sendCol = i_typecols(senderTypes, G, opts.condition);
recvCol = i_typecols(receiverTypes, G, opts.condition);

% -------------------------------------------------------------------------
% Apply EVERY matching Channel-B rule, multiplicatively
% -------------------------------------------------------------------------
% phi_ij = prod_k factor(f_k, a_k) over the modules k whose rule matches the
% pair, with a_k = +1 for an enabling glycan and -1 for a masking one. A pair
% may therefore carry several modules at once - sialyl-Lewis-x needs both
% sialylation and fucosylation - and that is what makes a per-pair dependency
% fingerprint expressible at all.
%
% This replaced first-matching-rule-wins, under which every rule after the
% first was silently dropped. With lambda = 0 each factor is exactly 1 and the
% product is 1, so the control identity is preserved regardless of how many
% rules matched.
detail = struct("pair", [], "module", strings(0, 1), ...
    "compartment", strings(0, 1), "exponent", [], "factor", [], ...
    "nodeGene", strings(0, 1), "nodeSites", [], "nodeFactor", []);

for j = 1:height(mapB)
    moduleRow = find(G.setnames == mapB.GlycoModule(j), 1);
    if isempty(moduleRow)
        continue;       % module not scored in this dataset; contributes 1
    end

    % cellstr(), not the string array itself: regexpi returns a bare double
    % (not wrapped in a cell) when its first argument has exactly one
    % element, so cellfun below errors on a single-pair call. cellstr()
    % forces cell-array semantics regardless of how many pairs there are.
    ligHit = ~cellfun(@isempty, regexpi(cellstr(ligNames), mapB.LigandPattern(j), "once"));
    recHit = ~cellfun(@isempty, regexpi(cellstr(recNames), mapB.ReceptorPattern(j), "once"));
    hit = ligHit(:) & recHit(:);
    if ~any(hit)
        continue;
    end

    if mapB.Compartment(j) == "sender"
        col = sendCol;
    else
        col = recvCol;
    end
    hit = hit & col > 0;        % skip cell types absent from G
    if ~any(hit)
        continue;
    end

    f = G.scores(moduleRow, col(hit))';
    a = mapB.Exponent(j);
    switch opts.mode
        case "center"
            factor = 1 + opts.lambda*a*(2*f - 1);
        case "boost"
            factor = 1 + opts.lambda*a*f;
        otherwise   % "raw": a masking module attenuates as capacity rises
            if a >= 0
                factor = f;
            else
                factor = 1 - f;
            end
    end

    idx = find(hit);

    % -- Node-level factor: breaks the tie within a symbol-pattern family ---
    % Same gene-selection convention as the module score above: sender
    % compartment reads the ligand (built by the sender cell), receiver
    % compartment reads the receptor (built by the receiver cell).
    if opts.usenode
        if mapB.Compartment(j) == "sender"
            nodeGene = ligNames(idx);
        else
            nodeGene = recNames(idx);
        end
        sites = i_nodesites(nodeGene, nodeT);
        [fNode, hasRange] = i_rangenorm(sites);
        if hasRange
            switch opts.mode
                case "center"
                    factorNode = 1 + opts.lambda*a*(2*fNode - 1);
                case "boost"
                    factorNode = 1 + opts.lambda*a*fNode;
                otherwise
                    if a >= 0
                        factorNode = fNode;
                    else
                        factorNode = 1 - fNode;
                    end
            end
        else
            % No distinct genes to discriminate among (a single gene matched,
            % or every matched gene ties on site count): factorNode is exactly
            % 1 regardless of mode. 0.5 is only "no effect" for "center"'s own
            % formula - "boost" (neutral at f=0) and "raw" (no neutral point at
            % all) would otherwise apply a real, unwarranted adjustment here.
            factorNode = ones(size(fNode));
        end
    else
        nodeGene = repmat("", numel(idx), 1);
        sites = nan(numel(idx), 1);
        factorNode = ones(numel(idx), 1);
    end

    w(hit) = w(hit).*factor.*factorNode;

    tag = mapB.GlycoModule(j);
    if a < 0
        tag = "-" + tag;
    end
    moduleName(idx) = i_joinstr(moduleName(idx), tag);
    compartment(idx) = i_joinstr(compartment(idx), mapB.Compartment(j));

    detail.pair = [detail.pair; idx];
    detail.module = [detail.module; repmat(mapB.GlycoModule(j), numel(idx), 1)];
    detail.compartment = [detail.compartment; repmat(mapB.Compartment(j), numel(idx), 1)];
    detail.exponent = [detail.exponent; repmat(a, numel(idx), 1)];
    detail.factor = [detail.factor; factor];
    detail.nodeGene = [detail.nodeGene; nodeGene(:)];
    detail.nodeSites = [detail.nodeSites; sites(:)];
    detail.nodeFactor = [detail.nodeFactor; factorNode(:)];
end

w = max(w, opts.minWeight);

end


%% ---- append a tag to each element of a string column, "+"-separated ----
function s = i_joinstr(s, tag)
empty = strlength(s) == 0;
s(empty) = tag;
s(~empty) = s(~empty) + "+" + tag;
end


%% ---- broadcast a scalar label to a full-length column ----
function v = i_broadcast(v, n, name)
if isscalar(v) && n > 1
    v = repmat(v, n, 1);
elseif numel(v) ~= n
    error("GLY:LRWEIGHT:BadSizes", ...
        "%s must be scalar or have %d elements (has %d).", name, n, numel(v));
end
end


%% ---- map each element's cell type to its column in G.scores ----
function col = i_typecols(types, G, condition)
[uTypes, ~, back] = unique(types);
uCol = zeros(numel(uTypes), 1);
for k = 1:numel(uTypes)
    uCol(k) = gly.groupcol(G, uTypes(k), condition);
end
col = uCol(back);
end


%% ---- total documented glycosite count per gene, 0 if not found ----
function sites = i_nodesites(genes, nodeT)
[tf, loc] = ismember(upper(genes(:)), nodeT.gene);
sites = zeros(numel(genes), 1);
sites(tf) = nodeT.n_total(loc(tf));
end


%% ---- min-max normalize to [0, 1]; flags whether there was any range ----
function [f, hasRange] = i_rangenorm(v)
v = v(:);
lo = min(v);
hi = max(v);
hasRange = hi > lo;
if hasRange
    f = (v - lo)./(hi - lo);
else
    f = 0.5*ones(size(v));   % placeholder; caller must not use it when ~hasRange
end
end
