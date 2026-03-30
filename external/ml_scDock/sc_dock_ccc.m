function result = sc_dock_ccc(X_norm, g, c_cell_type, varargin)
% SC_DOCK_CCC  Module 2 of scDock: cell-cell communication inference.
%
% Computes ligand-receptor (LR) communication probabilities between cell-type
% pairs using a permutation-based framework (analogous to CellChat v2).
% Returns ranked LR interactions for single-group or multi-group analysis.
%
% USAGE:
%   % Single-group
%   result = sc_dock_ccc(X_norm, g, c_cell_type)
%   result = sc_dock_ccc(X_norm, g, c_cell_type, 'lr_db', T_lr)
%
%   % Multi-group (compare two conditions)
%   result = sc_dock_ccc(X_norm, g, c_cell_type, ...
%               'condition', condition_vec, ...
%               'cond_labels', {'ctrl','treat'})
%
% INPUTS:
%   X_norm      - genes x cells normalized log1p expression matrix
%   g           - gene list (string array), length = size(X_norm,1)
%   c_cell_type - cell-type label per cell (string array or cellstr)
%
% OPTIONAL NAME-VALUE PAIRS:
%   'lr_db'         - LR database table with variables 'ligand' and 'receptor'
%                     (string columns). If empty, loads default human CellChatDB
%                     from lr_db_human.mat in the same folder. (default [])
%   'condition'     - condition label per cell (string/cellstr). When provided,
%                     multi-group mode is activated. (default [])
%   'cond_labels'   - 1x2 cell of condition names to compare, e.g. {'ctrl','treat'}
%                     (required for multi-group; default {})
%   'n_perm'        - number of permutations for significance test (default 100)
%   'pval_cutoff'   - p-value threshold (default 0.05)
%   'min_cells'     - minimum cells per cell-type to include (default 10)
%   'sender'        - restrict sender cell types (string array; [] = all)
%   'receiver'      - restrict receiver cell types (string array; [] = all)
%
% OUTPUT:
%   result  - struct with fields:
%     .T_interactions - table of LR interactions, columns:
%                       ligand, receptor, sender, receiver,
%                       prob, pval  [, prob_cond1, prob_cond2, delta_prob]
%     .T_pathways     - table of pathway-level summary (if 'pathway' column
%                       present in lr_db): pathway, sender, receiver, prob_sum
%     .cell_types     - unique cell-type labels used
%     .mode           - 'single' or 'multi'
%
% DEPENDENCIES (scGEAToolbox_dev):
%   pkg.e_fdr_bh, pkg.i_grpmean (via i_grpmean local fallback)
%
% See also: SC_DOCK_RNA, SC_DOCK_VINA, SC_DOCK

p = inputParser;
addRequired(p, 'X_norm',      @isnumeric);
addRequired(p, 'g');
addRequired(p, 'c_cell_type');
addParameter(p, 'lr_db',       [],     @(x) isempty(x) || istable(x));
addParameter(p, 'condition',   [],     @(x) isempty(x) || isvector(x));
addParameter(p, 'cond_labels', {},     @iscell);
addParameter(p, 'n_perm',      100,    @(x) isscalar(x) && x >= 1);
addParameter(p, 'pval_cutoff', 0.05,   @(x) isscalar(x) && x > 0 && x <= 1);
addParameter(p, 'min_cells',   10,     @(x) isscalar(x) && x >= 1);
addParameter(p, 'sender',      [],     @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(p, 'receiver',    [],     @(x) isempty(x) || isstring(x) || iscellstr(x));
parse(p, X_norm, g, c_cell_type, varargin{:});
opt = p.Results;

if iscellstr(g),           g           = string(g);           end %#ok<ISCLSTR>
if iscellstr(c_cell_type), c_cell_type = string(c_cell_type); end %#ok<ISCLSTR>

% -------------------------------------------------------------------------
% Load LR database
% -------------------------------------------------------------------------
if isempty(opt.lr_db)
    T_lr = i_load_default_lr_db();
else
    T_lr = opt.lr_db;
end
if ~all(ismember({'ligand','receptor'}, T_lr.Properties.VariableNames))
    error('sc_dock_ccc:BadLRDB', ...
        'lr_db table must have columns named ''ligand'' and ''receptor''.');
end
T_lr.ligand   = string(T_lr.ligand);
T_lr.receptor = string(T_lr.receptor);

% -------------------------------------------------------------------------
% Cell-type setup
% -------------------------------------------------------------------------
cell_types = unique(c_cell_type, 'stable');
% Filter by min_cells
keep_ct = false(numel(cell_types), 1);
for ci = 1:numel(cell_types)
    keep_ct(ci) = sum(c_cell_type == cell_types(ci)) >= opt.min_cells;
end
cell_types = cell_types(keep_ct);
if ~isempty(opt.sender),   sender_cts   = string(opt.sender);   else, sender_cts   = cell_types; end
if ~isempty(opt.receiver), receiver_cts = string(opt.receiver); else, receiver_cts = cell_types; end

fprintf('[sc_dock_ccc] %d cell types, %d LR pairs in database.\n', ...
    numel(cell_types), height(T_lr));

% -------------------------------------------------------------------------
% Determine mode
% -------------------------------------------------------------------------
multi_mode = ~isempty(opt.condition);
if multi_mode
    if numel(opt.cond_labels) ~= 2
        error('sc_dock_ccc:BadCondLabels', ...
            '''cond_labels'' must be a 1x2 cell, e.g. {''ctrl'',''treat''}.');
    end
    cond = string(opt.condition);
    cond1 = string(opt.cond_labels{1});
    cond2 = string(opt.cond_labels{2});
end

% -------------------------------------------------------------------------
% Pre-compute per-cell-type mean expression and cell index lookup
% -------------------------------------------------------------------------
fprintf('[sc_dock_ccc] Computing cell-type mean expression...\n');
n_ct = numel(cell_types);
ct_mean = zeros(numel(g), n_ct);
ct_cells = cell(n_ct, 1);   % cell indices per type (for permutation)
for ci = 1:n_ct
    idx = find(c_cell_type == cell_types(ci));
    ct_cells{ci} = idx;
    ct_mean(:, ci) = mean(X_norm(:, idx), 2);
end

% Map cell-type name -> ct_mean column index
ct2col = containers.Map(cellstr(cell_types), num2cell(1:n_ct));

% Map gene names to row indices (for fast lookup)
[~, g_idx_map_keys, g_idx_map_vals] = unique(g, 'stable');
gene2row = containers.Map(cellstr(g(g_idx_map_keys)), g_idx_map_vals);

% -------------------------------------------------------------------------
% Filter LR pairs to genes present in data
% -------------------------------------------------------------------------
has_lig = isKey(gene2row, cellstr(T_lr.ligand));
has_rec = isKey(gene2row, cellstr(T_lr.receptor));
T_lr = T_lr(has_lig & has_rec, :);
fprintf('[sc_dock_ccc] %d LR pairs after filtering to expressed genes.\n', height(T_lr));

% -------------------------------------------------------------------------
% Main loop: compute communication probability + permutation p-values
% -------------------------------------------------------------------------
n_lr      = height(T_lr);
n_sender  = numel(sender_cts);
n_recv    = numel(receiver_cts);
n_records = n_lr * n_sender * n_recv;

lig_out   = strings(n_records, 1);
rec_out   = strings(n_records, 1);
send_out  = strings(n_records, 1);
recv_out  = strings(n_records, 1);
prob_out  = zeros(n_records, 1);
pval_out  = ones(n_records, 1);
if multi_mode
    prob_c1   = zeros(n_records, 1);
    prob_c2   = zeros(n_records, 1);
end

row = 0;
fprintf('[sc_dock_ccc] Computing communication probabilities (n_perm=%d)...\n', ...
    opt.n_perm);

for li = 1:n_lr
    lig_gene = T_lr.ligand(li);
    rec_gene = T_lr.receptor(li);
    lig_row  = gene2row(char(lig_gene));
    rec_row  = gene2row(char(rec_gene));

    for si = 1:n_sender
        s_ct  = sender_cts(si);
        if ~isKey(ct2col, char(s_ct)), continue; end
        s_col  = ct2col(char(s_ct));
        s_cells = ct_cells{s_col};
        if numel(s_cells) < opt.min_cells, continue; end

        for ri = 1:n_recv
            r_ct  = receiver_cts(ri);
            if ~isKey(ct2col, char(r_ct)), continue; end
            r_col  = ct2col(char(r_ct));
            r_cells = ct_cells{r_col};
            if numel(r_cells) < opt.min_cells, continue; end

            row = row + 1;
            lig_out(row)  = lig_gene;
            rec_out(row)  = rec_gene;
            send_out(row) = s_ct;
            recv_out(row) = r_ct;

            % Observed probability from pre-computed means (no recomputation)
            lig_mean_s = ct_mean(lig_row, s_col);
            rec_mean_r = ct_mean(rec_row, r_col);
            prob = lig_mean_s * rec_mean_r;

            if multi_mode
                [prob, pval, p1, p2] = i_comm_prob_multi( ...
                    X_norm, lig_row, rec_row, s_cells, r_cells, ...
                    cond, cond1, cond2, opt.n_perm);
                prob_c1(row) = p1;
                prob_c2(row) = p2;
            else
                pval = i_perm_pval(X_norm, lig_row, rec_row, ...
                    s_cells, r_cells, prob, opt.n_perm);
            end
            prob_out(row) = prob;
            pval_out(row) = pval;
        end
    end
end

% Trim to actual rows filled
lig_out  = lig_out(1:row);
rec_out  = rec_out(1:row);
send_out = send_out(1:row);
recv_out = recv_out(1:row);
prob_out = prob_out(1:row);
pval_out = pval_out(1:row);

% -------------------------------------------------------------------------
% FDR correction and filtering
% -------------------------------------------------------------------------
[~, ~, ~, adj_p] = pkg.e_fdr_bh(pval_out, opt.pval_cutoff, 'pdep', 'no');
sig_idx = adj_p <= opt.pval_cutoff;

if multi_mode
    T_int = table(lig_out, rec_out, send_out, recv_out, prob_out, pval_out, adj_p, ...
        prob_c1(1:row), prob_c2(1:row), ...
        prob_c2(1:row) - prob_c1(1:row), ...
        'VariableNames', {'ligand','receptor','sender','receiver', ...
                          'prob','pval','adj_pval', ...
                          'prob_cond1','prob_cond2','delta_prob'});
    T_int = T_int(sig_idx, :);
    % Rank by |delta_prob| descending
    T_int = sortrows(T_int, 'delta_prob', 'descend');
else
    T_int = table(lig_out, rec_out, send_out, recv_out, prob_out, pval_out, adj_p, ...
        'VariableNames', {'ligand','receptor','sender','receiver', ...
                          'prob','pval','adj_pval'});
    T_int = T_int(sig_idx, :);
    % Rank by prob descending
    T_int = sortrows(T_int, 'prob', 'descend');
end

fprintf('[sc_dock_ccc] %d significant LR interactions (FDR < %.2f).\n', ...
    height(T_int), opt.pval_cutoff);

% -------------------------------------------------------------------------
% Pathway-level summary (if 'pathway' column present)
% -------------------------------------------------------------------------
T_path = table();
if ismember('pathway', T_lr.Properties.VariableNames) && height(T_int) > 0
    T_lr.pathway = string(T_lr.pathway);
    % Join pathway info onto T_int
    [~, lr_row] = ismember( ...
        strcat(T_int.ligand, '_', T_int.receptor), ...
        strcat(T_lr.ligand,  '_', T_lr.receptor));
    pathway_col = repmat("", height(T_int), 1);
    has_match = lr_row > 0;
    pathway_col(has_match) = T_lr.pathway(lr_row(has_match));
    T_int.pathway = pathway_col;

    % Aggregate
    [pw_keys, ~, pw_ic] = unique( ...
        strcat(pathway_col, '|', T_int.sender, '|', T_int.receiver));
    pw_prob = accumarray(pw_ic, T_int.prob, [], @sum);
    parts = split(pw_keys, '|');
    T_path = table(parts(:,1), parts(:,2), parts(:,3), pw_prob, ...
        'VariableNames', {'pathway','sender','receiver','prob_sum'});
    T_path = sortrows(T_path, 'prob_sum', 'descend');
end

% -------------------------------------------------------------------------
% Pack results
% -------------------------------------------------------------------------
result.T_interactions = T_int;
result.T_pathways     = T_path;
result.cell_types     = cell_types;
result.mode           = 'single';
if multi_mode, result.mode = 'multi'; end

fprintf('[sc_dock_ccc] Done.\n');
end

% =========================================================================
% Local helpers
% =========================================================================

function T_lr = i_load_default_lr_db()
% Load LR database from scGEAToolbox_dev assets (searched on MATLAB path
% and relative to this file). Prefers Ligand_Receptor_more.txt (cleanest
% format), falls back to Ligand_Receptor.txt, then connectomeDB2020.txt.
asset_dir = i_find_lr_asset_dir();

% Option 0: CellChatDB export (richest — includes pathway column)
f0 = fullfile(asset_dir, 'CellChatDB_human_LR.txt');
if isfile(f0)
    T_lr = readtable(f0, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    return
end

% Option 1: Ligand_Receptor_more.txt  (columns: ligand, receptor)
f1 = fullfile(asset_dir, 'Ligand_Receptor_more.txt');
if isfile(f1)
    T_lr = readtable(f1, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    T_lr = T_lr(:, {'ligand','receptor'});
    return
end

% Option 2: Ligand_Receptor.txt  (columns include Ligand.ApprovedSymbol, Receptor.ApprovedSymbol)
f2 = fullfile(asset_dir, 'Ligand_Receptor.txt');
if isfile(f2)
    T_raw = readtable(f2, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    T_lr = table(T_raw.("Ligand_ApprovedSymbol"), T_raw.("Receptor_ApprovedSymbol"), ...
        'VariableNames', {'ligand','receptor'});
    return
end

% Option 3: connectomeDB2020.txt  (columns: Ligand gene symbol, Receptor gene symbol)
f3 = fullfile(asset_dir, 'connectomeDB2020.txt');
if isfile(f3)
    T_raw = readtable(f3, 'FileType', 'text', 'Delimiter', '\t', ...
        'TextType', 'string');
    T_lr = table(T_raw.("Ligand_gene_symbol"), T_raw.("Receptor_gene_symbol"), ...
        'VariableNames', {'ligand','receptor'});
    return
end

% Option 4: .mat files
f4 = fullfile(asset_dir, 'Ligand_Receptor_more.mat');
if isfile(f4)
    S = load(f4);
    fn = fieldnames(S);
    T_lr = S.(fn{1});
    if ~istable(T_lr)
        % Assume Nx2 char/string array
        T_lr = array2table(T_lr, 'VariableNames', {'ligand','receptor'});
    end
    return
end

error('sc_dock_ccc:NoLRDB', ...
    ['Cannot find LR database. Ensure scGEAToolbox_dev is on the MATLAB path, ' ...
     'or supply a table via the ''lr_db'' parameter with columns ' ...
     '''ligand'' and ''receptor''.']);
end

function asset_dir = i_find_lr_asset_dir()
% Search for the Ligand_Receptor asset folder
candidates = {
    fullfile(fileparts(mfilename('fullpath')), ...
        '..', 'scGEAToolbox_dev', 'assets', 'Ligand_Receptor');
    fullfile(fileparts(mfilename('fullpath')), ...
        '..', '..', 'scGEAToolbox_dev', 'assets', 'Ligand_Receptor');
};
% Also search MATLAB path for the asset folder
path_dirs = strsplit(path, pathsep);
for k = 1:numel(path_dirs)
    if contains(path_dirs{k}, 'scGEAToolbox')
        candidates{end+1} = fullfile(fileparts(path_dirs{k}), ...
            'assets', 'Ligand_Receptor'); %#ok<AGROW>
        candidates{end+1} = fullfile(path_dirs{k}, ...
            'assets', 'Ligand_Receptor'); %#ok<AGROW>
    end
end
for k = 1:numel(candidates)
    d = candidates{k};
    if isfolder(d) && isfile(fullfile(d, 'Ligand_Receptor_more.txt'))
        asset_dir = d;
        return
    end
end
% Last resort: hardcoded sibling path
asset_dir = 'C:\Users\jingc\Documents\GitHub\scGEAToolbox_dev\assets\Ligand_Receptor';
end

function pval = i_perm_pval(X, lig_row, rec_row, s_cells, r_cells, obs_prob, n_perm)
% Permutation p-value for single-group communication probability.
% obs_prob already computed from ct_mean; only null distribution needed here.
% Shuffles cell-type membership by drawing random same-size subsets.
n_cells = size(X, 2);
n_s = numel(s_cells);
n_r = numel(r_cells);
null_probs = zeros(n_perm, 1);
for b = 1:n_perm
    perm    = randperm(n_cells);
    s_p     = perm(1:n_s);
    r_p     = perm(n_s+1:n_s+n_r);
    null_probs(b) = mean(X(lig_row, s_p)) * mean(X(rec_row, r_p));
end
pval = mean(null_probs >= obs_prob);
end

function [prob, pval, p1, p2] = i_comm_prob_multi( ...
        X, lig_row, rec_row, s_cells, r_cells, ...
        cond, cond1, cond2, n_perm)
% Compute per-condition probabilities and permutation p-value on |delta|.
% s_cells / r_cells are index vectors (not logical masks).
s1 = s_cells(cond(s_cells) == cond1);
r1 = r_cells(cond(r_cells) == cond1);
s2 = s_cells(cond(s_cells) == cond2);
r2 = r_cells(cond(r_cells) == cond2);

p1 = 0; p2 = 0;
if ~isempty(s1) && ~isempty(r1)
    p1 = mean(X(lig_row, s1)) * mean(X(rec_row, r1));
end
if ~isempty(s2) && ~isempty(r2)
    p2 = mean(X(lig_row, s2)) * mean(X(rec_row, r2));
end
prob = abs(p2 - p1);

% Permutation: shuffle condition labels within each cell-type pool
n_s = numel(s_cells);  n_r = numel(r_cells);
n_s1 = numel(s1);      n_r1 = numel(r1);

null_delta = zeros(n_perm, 1);
for b = 1:n_perm
    ps = randperm(n_s);
    pr = randperm(n_r);
    sp1 = s_cells(ps(1:n_s1));
    rp1 = r_cells(pr(1:n_r1));
    sp2 = s_cells(ps(n_s1+1:end));
    rp2 = r_cells(pr(n_r1+1:end));
    np1 = 0; np2 = 0;
    if ~isempty(sp1) && ~isempty(rp1)
        np1 = mean(X(lig_row, sp1)) * mean(X(rec_row, rp1));
    end
    if ~isempty(sp2) && ~isempty(rp2)
        np2 = mean(X(lig_row, sp2)) * mean(X(rec_row, rp2));
    end
    null_delta(b) = abs(np2 - np1);
end
pval = mean(null_delta >= prob);
end
