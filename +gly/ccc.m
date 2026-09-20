function result = ccc(X, genelist, c_cell_type, opts)
%GLY.CCC  Glyco-lectin cell-cell communication.
%   result = GLY.CCC(X, genelist, c_cell_type) infers glycan-mediated
%   communication between cell types. Unlike protein ligand-receptor engines,
%   the "ligand" here is a glycan determinant that is invisible to RNA-seq: its
%   abundance is approximated by the sender's biosynthetic capacity (a per-cell
%   glyco-module score from GLY.STATE), and the "receptor" is a lectin gene
%   whose expression is measured on the receiver. Each sender->receiver edge is
%   scored as (mean sender glyco-module score) x (mean receiver lectin
%   expression) and tested against a permutation null, mirroring SC_DOCK_CCC.
%
%   The glycan-to-lectin cognate map comes from GLY.LECTINMAP and the
%   lectin gene memberships from GLY.GENESETS.
%
%   USAGE:
%     result = gly.ccc(sce.X, sce.g, sce.c_cell_type_tx);
%     result = gly.ccc(X, genelist, c_cell_type, n_perm=200);
%
%   INPUTS:
%     X           - genes-by-cells expression matrix (normalized, e.g. log1p)
%     genelist    - G-by-1 gene symbols (length = rows of X)
%     c_cell_type - N-by-1 cell-type label per cell (length = columns of X)
%
%   NAME-VALUE ARGUMENTS:
%     methodid    - scoring method for GLY.STATE (1 UCell, 2 AddModuleScore
%                   default, 3 AUCell)
%     minGenes    - minimum module genes present to score a module (default 3)
%     min_cells   - minimum cells per cell type to include (default 10)
%     n_perm      - number of permutations (default 100)
%     pval_cutoff - FDR threshold for the returned table (default 0.05)
%     sender      - restrict sender cell types (string array; [] = all)
%     receiver    - restrict receiver cell types (string array; [] = all)
%
%   OUTPUT:
%     result - struct with fields:
%       .T_interactions - table of significant edges, columns:
%                         glyco_module, lectin, sender, receiver,
%                         prob, pval, adj_pval (ranked by prob descending)
%       .cell_types     - cell-type labels used
%       .epitopes       - Channel-A epitope rows that contributed
%
%   Because per-cell glyco-module scores are relative (centred on background),
%   negative scores are clipped to zero before use as a ligand signal (below-
%   background biosynthesis contributes no glycan ligand).
%
% see also: GLY.STATE, GLY.LECTINMAP, GLY.GENESETS,
%           GLY.WEIGHT, SC_DOCK_CCC

arguments
    X {mustBeNumeric}
    genelist (:, 1) string
    c_cell_type (:, 1) string
    opts.methodid (1, 1) double = 2
    opts.minGenes (1, 1) double = 3
    opts.min_cells (1, 1) double = 10
    opts.n_perm (1, 1) double = 100
    opts.pval_cutoff (1, 1) double = 0.05
    opts.sender (:, 1) string = strings(0, 1)
    opts.receiver (:, 1) string = strings(0, 1)
end

if numel(genelist) ~= size(X, 1)
    error("GLY:CCC:BadGenelist", ...
        "GENELIST length (%d) must equal the number of rows of X (%d).", ...
        numel(genelist), size(X, 1));
end
if numel(c_cell_type) ~= size(X, 2)
    error("GLY:CCC:BadLabels", ...
        "C_CELL_TYPE length (%d) must equal the number of columns of X (%d).", ...
        numel(c_cell_type), size(X, 2));
end

% -------------------------------------------------------------------------
% Per-cell glyco-module scores (sender ligand signal) and lectin memberships
% -------------------------------------------------------------------------
[cs, setnames] = gly.state(X, genelist, opts.methodid, opts.minGenes);
[setmatrx, allnames, setgenes] = gly.genesets();
mapA = gly.lectinmap();

% -------------------------------------------------------------------------
% Cell-type setup (filter by min_cells)
% -------------------------------------------------------------------------
cell_types = unique(c_cell_type, "stable");
keep_ct = false(numel(cell_types), 1);
for ci = 1:numel(cell_types)
    keep_ct(ci) = sum(c_cell_type == cell_types(ci)) >= opts.min_cells;
end
cell_types = cell_types(keep_ct);
if isempty(cell_types)
    error("GLY:CCC:NoCellTypes", ...
        "No cell type has at least min_cells (%d) cells.", opts.min_cells);
end

if isempty(opts.sender),   sender_cts   = cell_types; else, sender_cts   = opts.sender;   end
if isempty(opts.receiver), receiver_cts = cell_types; else, receiver_cts = opts.receiver; end
sender_cts   = sender_cts(ismember(sender_cts, cell_types));
receiver_cts = receiver_cts(ismember(receiver_cts, cell_types));

ct_cells = cell(numel(cell_types), 1);
for ci = 1:numel(cell_types)
    ct_cells{ci} = find(c_cell_type == cell_types(ci));
end
ct2idx = containers.Map(cellstr(cell_types), num2cell(1:numel(cell_types)));

% -------------------------------------------------------------------------
% Pre-pass: build the ligand signal and present lectin genes per epitope
% -------------------------------------------------------------------------
geneUpper = upper(genelist);
nEpi = height(mapA);
ligsigCell   = cell(nEpi, 1);       % 1-by-N clipped ligand signal per epitope
modLabelStr  = strings(nEpi, 1);    % combined module label per epitope
lecRowsCell  = cell(nEpi, 1);       % gene-row indices of present lectins
lecNameCell  = cell(nEpi, 1);       % gene names of present lectins
epiValid     = false(nEpi, 1);

for ei = 1:nEpi
    reqModules = pkg.i_str2genelist(mapA.GlycoModules(ei));
    [inScore, scoreRows] = ismember(reqModules, setnames);
    if ~any(inScore), continue; end

    lectinModule = mapA.LectinModule(ei);
    lmRow = strcmp(allnames, lectinModule);
    if ~any(lmRow), continue; end
    lectinGenes = setgenes(setmatrx(lmRow, :));
    [hasLec, lecRows] = ismember(upper(lectinGenes), geneUpper);
    if ~any(hasLec), continue; end

    ligsig = mean(cs(scoreRows(inScore), :), 1);
    ligsigCell{ei}  = max(ligsig, 0);       % clip below-background to zero
    modLabelStr(ei) = strjoin(reqModules(inScore), "+");
    lecRowsCell{ei} = lecRows(hasLec);
    lecNameCell{ei} = lectinGenes(hasLec);
    epiValid(ei)    = true;
end

% -------------------------------------------------------------------------
% Pre-allocate to the exact number of edges, then fill
% -------------------------------------------------------------------------
nPairs = numel(sender_cts) * numel(receiver_cts);
nEdges = 0;
for ei = 1:nEpi
    if epiValid(ei)
        nEdges = nEdges + numel(lecRowsCell{ei}) * nPairs;
    end
end

gm_out   = strings(nEdges, 1);
lec_out  = strings(nEdges, 1);
send_out = strings(nEdges, 1);
recv_out = strings(nEdges, 1);
prob_out = zeros(nEdges, 1);
pval_out = zeros(nEdges, 1);
epitopes_used = mapA.Epitope(epiValid);

row = 0;
for ei = 1:nEpi
    if ~epiValid(ei), continue; end
    ligsig   = ligsigCell{ei};
    modLabel = modLabelStr(ei);
    lecRows  = lecRowsCell{ei};
    lecNames = lecNameCell{ei};

    for li = 1:numel(lecRows)
        recvec = full(X(lecRows(li), :));
        for si = 1:numel(sender_cts)
            s_cells = ct_cells{ct2idx(char(sender_cts(si)))};
            lig_mean_s = mean(ligsig(s_cells));
            for ri = 1:numel(receiver_cts)
                r_cells = ct_cells{ct2idx(char(receiver_cts(ri)))};
                obs = lig_mean_s * mean(recvec(r_cells));
                pval = i_perm_pval(ligsig, recvec, s_cells, r_cells, ...
                    obs, opts.n_perm);

                row = row + 1;
                gm_out(row)   = modLabel;
                lec_out(row)  = lecNames(li);
                send_out(row) = sender_cts(si);
                recv_out(row) = receiver_cts(ri);
                prob_out(row) = obs;
                pval_out(row) = pval;
            end
        end
    end
end

% -------------------------------------------------------------------------
% FDR correction, filtering, ranking
% -------------------------------------------------------------------------
if row == 0
    warning("GLY:CCC:NoEdges", ...
        "No glyco-lectin edges could be formed from the data.");
    T_int = table(Size=[0 7], ...
        VariableTypes=["string", "string", "string", "string", ...
        "double", "double", "double"], ...
        VariableNames=["glyco_module", "lectin", "sender", "receiver", ...
        "prob", "pval", "adj_pval"]);
else
    adj_p = pkg.e_fdr(pval_out);
    T_int = table(gm_out, lec_out, send_out, recv_out, prob_out, pval_out, ...
        adj_p, VariableNames=["glyco_module", "lectin", "sender", ...
        "receiver", "prob", "pval", "adj_pval"]);
    T_int = T_int(adj_p <= opts.pval_cutoff, :);
    T_int = sortrows(T_int, "prob", "descend");
end

result.T_interactions = T_int;
result.cell_types     = cell_types;
result.epitopes       = epitopes_used;

end


%% ---- permutation p-value for a glyco-ligand x lectin-receptor edge ----
function pval = i_perm_pval(ligvec, recvec, s_cells, r_cells, obs, n_perm)
% Null: draw random same-size sender and receiver groups from all cells and
% recompute the product of mean ligand signal and mean receptor expression.
% The two groups are drawn independently so the null is valid even for
% autocrine edges (sender and receiver the same cell type).
n_cells = numel(ligvec);
n_s = numel(s_cells);
n_r = numel(r_cells);
count = 0;
for b = 1:n_perm
    s_p = randperm(n_cells, n_s);
    r_p = randperm(n_cells, n_r);
    null_prob = mean(ligvec(s_p)) * mean(recvec(r_p));
    if null_prob >= obs
        count = count + 1;
    end
end
% (count + 1)/(n_perm + 1), not count/n_perm. The plain ratio is zero
% whenever no permutation beats the observation, and a zero survives BH as
% zero -- PKG.E_FDR([0 0.01 0.2 0.5]) returns [0 0.02 0.267 0.5] -- so such
% an edge passes any cutoff whatever the family size, while an edge one
% permutation worse (1/n_perm) is filtered out. With the default n_perm =
% 100 that made membership of the returned table a coin flip. TEN.TNGRN
% already uses the +1 form for the same reason.
pval = (count + 1) / (n_perm + 1);
end
