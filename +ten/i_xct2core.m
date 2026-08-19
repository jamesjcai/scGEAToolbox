function T = i_xct2core(sce1, sce2, celltype1, celltype2, twosided, cfg)
%I_XCT2CORE  Shared engine for the differential scTenifoldXct entry points.
%   T = TEN.I_XCT2CORE(SCE1, SCE2, CELLTYPE1, CELLTYPE2, TWOSIDED, CFG) runs the
%   two-sample comparison: four GRNs, one correspondence block per sample, a
%   spectral embedding per sample, and a chi-square test on the squared
%   difference of the aligned distances.
%
%   This stands to TEN.SCTENIFOLDXCT2 as TEN.I_XCTCORE stands to
%   TEN.SCTENIFOLDXCT: everything common lives here, and the only thing an entry
%   point may vary is the content of the correspondence block, supplied through
%   cfg.corrfcn. This function has no knowledge of glycobiology.
%
%   INPUTS:
%     SCE1, SCE2  - SingleCellExperiment objects for the two conditions
%     CELLTYPE1   - source cell type label
%     CELLTYPE2   - target cell type label
%     TWOSIDED    - logical; when true both directions are run and T is {T1, T2}
%
%   CFG FIELDS:
%     n_dim      - spectral embedding dimensions
%     mu         - cross-cell-type correspondence weight
%     verbose    - print progress
%     w12mode    - "lr" (cognate pairs) or "outer" (reference outer product)
%     alpha      - mean/variance blend for w12mode="outer"
%     grnoffset  - constant added to every GRN element; see TEN.I_XCTBLOCK
%     useparallel- pass a parfor down to net.pcrnet; see TEN.I_XCTGRN
%     tag        - string used to prefix progress messages
%     corrfcn    - [] for the published weight-1 correspondence, or a handle
%                     [li, ri, w, module, channel] = ...
%                         corrfcn(g, li, ri, src, tgt, sampleIdx)
%                  called once per direction PER SAMPLE, with sampleIdx 1 or 2.
%     provenance - logical; when true the output gains channel, glyco_module,
%                  glyco_weight_s1 and glyco_weight_s2 columns.
%
%   OUTPUT:
%     T - table: ligand, receptor, dist_s1, dist_s2, diff2, FC, p_value,
%         q_value, sorted by diff2 descending; {T1, T2} when TWOSIDED.
%
%   THE TWO SAMPLES MUST SHARE ONE CANDIDATE SET. The test compares distances
%   measured in two separately-computed embeddings, which is only meaningful if
%   the same gene pairs are being read out of both. A corrfcn may therefore
%   reweight pairs per sample - that is the point of the sample index - but it
%   may not give the two samples different pairs. The index vectors it returns
%   are compared and a mismatch is an error rather than a silent truncation.
%
% see also: TEN.SCTENIFOLDXCT2, TEN.SCTENIFOLDXCT2_GLYCO, TEN.I_XCTCHI2,
%           TEN.I_XCTW12, TEN.I_XCTBLOCK

tag = cfg.tag;
verbose = cfg.verbose;

% -- Expression matrices --------------------------------------------------
[X_s1, X_t1, X_s2, X_t2, g] = i_prepare(sce1, sce2, celltype1, celltype2, ...
    tag, verbose);

% -- Ligand-receptor database ---------------------------------------------
[lig_db, rec_db] = i_loadlrdb();
if verbose
    fprintf('[%s] L-R database: %d pairs loaded.\n', tag, numel(lig_db));
end

[li_idx, ri_idx] = i_matchlr(upper(string(g(:))), lig_db, rec_db);
if isempty(li_idx)
    warning('ten:i_xct2core:noPairs', 'No L-R pairs found in gene list.');
    T = table();
    return;
end
if verbose
    fprintf('[%s] %d L-R pairs matched in data.\n', tag, numel(li_idx));
end

% -- GRNs -----------------------------------------------------------------
% Same shared builder and settings as ten.sctenifoldxct: ncomp=3 with a 0.75
% edge filter. TEN.I_XCTGRN owns the pcrnet call, the max-normalisation, the
% filtering and the symmetrisation.
if verbose, fprintf('[%s] Building GRNs for sample 1 ...\n', tag); end
A_s1 = ten.i_xctgrn(X_s1, 3, 0.75, false, cfg.useparallel);
A_t1 = ten.i_xctgrn(X_t1, 3, 0.75, false, cfg.useparallel);
if verbose, fprintf('[%s] Building GRNs for sample 2 ...\n', tag); end
A_s2 = ten.i_xctgrn(X_s2, 3, 0.75, false, cfg.useparallel);
A_t2 = ten.i_xctgrn(X_t2, 3, 0.75, false, cfg.useparallel);

% One bundle per sample, so the direction swap is a single operation and the
% argument list stays readable.
S1 = struct('X_s', X_s1, 'X_t', X_t1, 'A_s', A_s1, 'A_t', A_t1);
S2 = struct('X_s', X_s2, 'X_t', X_t2, 'A_s', A_s2, 'A_t', A_t2);

% -- Differential alignment, one direction at a time ----------------------
if verbose
    fprintf('[%s] Differential alignment %s -> %s ...\n', tag, ...
        celltype1, celltype2);
end
T1 = i_xct2dir(g, S1, S2, li_idx, ri_idx, ...
    string(celltype1), string(celltype2), cfg);

if twosided
    if verbose
        fprintf('[%s] Differential alignment %s -> %s ...\n', tag, ...
            celltype2, celltype1);
    end
    T2 = i_xct2dir(g, i_swapdir(S1), i_swapdir(S2), li_idx, ri_idx, ...
        string(celltype2), string(celltype1), cfg);
    T = {T1, T2};
else
    T = T1;
end

end % i_xct2core


%% ---- differential alignment for one direction (source -> target) ----
function Td = i_xct2dir(g, S1, S2, li_idx, ri_idx, src, tgt, cfg)

ng = size(S1.A_s, 1);
verbose = cfg.verbose;

% -- Correspondence weights, per sample -----------------------------------
% Default: every pair gets weight 1 in both samples, exactly as published.
n_valid = numel(li_idx);
li1 = li_idx;
ri1 = ri_idx;
w1 = ones(n_valid, 1);
w2 = ones(n_valid, 1);
lr_module = repmat("", n_valid, 1);
channel = repmat("protein_LR", n_valid, 1);

if ~isempty(cfg.corrfcn)
    [li1, ri1, w1, lr_module, channel] = ...
        cfg.corrfcn(g, li_idx, ri_idx, src, tgt, 1);
    [li2, ri2, w2, ~, ~] = cfg.corrfcn(g, li_idx, ri_idx, src, tgt, 2);
    if ~isequal(li1, li2) || ~isequal(ri1, ri2)
        error("TEN:I_XCT2CORE:CandidateSetMismatch", ...
            "The correspondence hook returned different gene pairs for the " + ...
            "two samples (%d and %d pairs). The differential test compares " + ...
            "distances read out of two separate embeddings, so both must be " + ...
            "read at the same pairs. Vary the weights, not the pairs.", ...
            numel(li1), numel(li2));
    end
end

% -- Embeddings, one per sample -------------------------------------------
[P_s1, P_t1] = i_embed(S1, li1, ri1, w1, ng, cfg);
[P_s2, P_t2] = i_embed(S2, li1, ri1, w2, ng, cfg);

% All-pairs distances in each sample. The chi-square test scales by the mean
% statistic over EVERY pair and applies FDR across all of them, exactly as
% merge.py builds its full table before calling chi2_diff_test, so the
% candidate-only distances are not sufficient here.
D1 = pdist2(double(P_s1), double(P_t1));
D2 = pdist2(double(P_s2), double(P_t2));

diff2_all = (D1 - D2).^2;                       % merge.py:110
[p_all, q_all, FC_all] = ten.i_xctchi2(diff2_all(:), dof=1, tail="right");

candLin = sub2ind([ng, ng], li1, ri1);
d1 = D1(candLin);
d2 = D2(candLin);

% Output table, restricted to candidates after the correction as stat.py does
lig_out = g(li1);
rec_out = g(ri1);
Td = table(lig_out, rec_out, d1, d2, diff2_all(candLin), FC_all(candLin), ...
    p_all(candLin), q_all(candLin), ...
    'VariableNames', {'ligand', 'receptor', 'dist_s1', 'dist_s2', 'diff2', ...
    'FC', 'p_value', 'q_value'});

% Provenance columns are appended only on request, so the default output
% schema is exactly the published one. Both weights are carried, because in the
% differential setting their DIFFERENCE is the part that can move diff2.
if cfg.provenance
    Td.channel = channel;
    Td.glyco_module = lr_module;
    Td.glyco_weight_s1 = w1;
    Td.glyco_weight_s2 = w2;
end

Td = sortrows(Td, 'diff2', 'descend');          % stat.py:89

if verbose
    fprintf('[%s]   %d L-R pairs; %d with q < 0.05.\n', ...
        cfg.tag, height(Td), sum(Td.q_value < 0.05));
end

end % i_xct2dir


%% ---- spectral manifold alignment for one sample ----
function [P_s, P_t] = i_embed(S, li_idx, ri_idx, w, ng, cfg)
% Correspondence via TEN.I_XCTW12 and assembly via TEN.I_XCTBLOCK, so this
% shares one implementation with the other entry points.
%
% In "outer" mode the correspondence is a dense expression outer product in
% which the L-R database plays no part, so W is ignored - the same convention
% TEN.I_XCTCORE follows.

if cfg.w12mode == "outer"
    W12 = ten.i_xctw12(S.X_s, S.X_t, ng, mode="outer", alpha=cfg.alpha);
else
    W12 = ten.i_xctw12(S.X_s, S.X_t, ng, mode="lr", ...
        li_idx=li_idx, ri_idx=ri_idx, weights=w);
end

W = ten.i_xctblock(0.5*(S.A_s + S.A_s'), 0.5*(S.A_t + S.A_t'), W12, ...
    cfg.mu, cfg.grnoffset);

[P_s, P_t] = ten.i_laplacianembed(W, ng, cfg.n_dim, false);

end % i_embed


%% ---- reverse the source/target roles within one sample ----
function S = i_swapdir(S)
[S.X_s, S.X_t] = deal(S.X_t, S.X_s);
[S.A_s, S.A_t] = deal(S.A_t, S.A_s);
end % i_swapdir


%% ---- expression matrices for both samples ----
function [X_s1, X_t1, X_s2, X_t2, g] = i_prepare(sce1, sce2, ct1, ct2, ...
    tag, verbose)
[X_s1, X_t1] = i_fromsce(sce1, ct1, ct2);
[X_s2, X_t2] = i_fromsce(sce2, ct1, ct2);
g = sce1.g;

if isempty(X_s1) || isempty(X_t1) || isempty(X_s2) || isempty(X_t2)
    error("TEN:I_XCT2CORE:MissingCellType", ...
        "Cell counts are %d/%d in sample 1 and %d/%d in sample 2 for " + ...
        """%s""/""%s"". Both cell types must be present in both samples.", ...
        size(X_s1, 2), size(X_t1, 2), size(X_s2, 2), size(X_t2, 2), ct1, ct2);
end
if numel(g) ~= numel(sce2.g) || ~all(string(g(:)) == string(sce2.g(:)))
    error("TEN:I_XCT2CORE:GeneMismatch", ...
        "SCE1 and SCE2 must carry the same gene list in the same order " + ...
        "(%d and %d genes). Subset both to their shared genes first.", ...
        numel(g), numel(sce2.g));
end

if verbose
    fprintf('[%s] S1  %s: %d cells  |  %s: %d cells\n', tag, ...
        ct1, size(X_s1, 2), ct2, size(X_t1, 2));
    fprintf('[%s] S2  %s: %d cells  |  %s: %d cells\n', tag, ...
        ct1, size(X_s2, 2), ct2, size(X_t2, 2));
end
end % i_prepare


function [Xs, Xt] = i_fromsce(sce, c1, c2)
X = sce.X;
if issparse(X), X = full(X); end
X = single(X);
Xs = X(:, sce.c_cell_type_tx == c1);
Xt = X(:, sce.c_cell_type_tx == c2);
end % i_fromsce


%% ---- indices of database pairs present in the gene list ----
function [li_idx, ri_idx] = i_matchlr(g_up, lig_db, rec_db)
n_lr = numel(lig_db);
li = zeros(n_lr, 1);
ri = zeros(n_lr, 1);
n = 0;
for k = 1:n_lr
    l = find(g_up == lig_db(k), 1);
    r = find(g_up == rec_db(k), 1);
    if ~isempty(l) && ~isempty(r)
        n = n + 1;
        li(n) = l;
        ri(n) = r;
    end
end
li_idx = li(1:n);
ri_idx = ri(1:n);
end % i_matchlr


%% ---- built-in ligand-receptor database ----
function [lig_db, rec_db] = i_loadlrdb()
pw1 = fileparts(mfilename('fullpath'));
lr_mat = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.mat');
lr_txt = fullfile(pw1, '..', 'assets', 'Ligand_Receptor', 'Ligand_Receptor.txt');

if exist(lr_mat, 'file')
    db = load(lr_mat, 'ligand', 'receptor');
    lig_db = upper(string(db.ligand(:)));
    rec_db = upper(string(db.receptor(:)));
else
    T_lr = readtable(lr_txt, 'FileType', 'text', 'Delimiter', '\t');
    lig_db = upper(string(T_lr{:, 3}));
    rec_db = upper(string(T_lr{:, 5}));
end
end % i_loadlrdb
