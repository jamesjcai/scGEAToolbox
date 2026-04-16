function result = sc_dock(X, g, varargin)
% SC_DOCK  scDock pipeline: scRNA-seq → cell-cell communication → docking.
%
% End-to-end reimplementation of scDock (Bioinformatics btag103) in MATLAB.
% Runs three sequential modules:
%   1. scRNA-seq preprocessing and cell-type annotation  (sc_dock_rna)
%   2. Cell-cell communication inference                  (sc_dock_ccc)
%   3. AutoDock Vina molecular docking                    (sc_dock_vina)
%
% USAGE:
%   % Minimal — only run scRNA-seq module
%   result = sc_dock(X, g)
%
%   % Full pipeline with compounds
%   result = sc_dock(X, g, ...
%       'compounds', {'aspirin','ibuprofen'}, ...
%       'species',   'human', ...
%       'outdir',    './my_dock_run')
%
%   % Skip docking (no compounds provided)
%   result = sc_dock(X, g, 'compounds', [])
%
% REQUIRED INPUTS:
%   X  - genes x cells raw count matrix (sparse or dense)
%   g  - gene list (string array or cell array of chars)
%
% OPTIONAL NAME-VALUE PAIRS (forwarded to sub-modules):
%
%  --- Module 1 (sc_dock_rna) ---
%   'batch_id'        - batch label vector (default [])
%   'species'         - 'human' (default) | 'mouse'
%   'libsz_cutoff'    - min library size (default 500)
%   'mt_ratio'        - max mt read ratio (default 0.15)
%   'min_cells'       - min cells per gene (default 10)
%   'n_hvg'           - number of HVGs (default 2000)
%   'n_pcs'           - PCs to use; 0 = auto elbow (default 0)
%   'n_clusters'      - clusters for sc_cluster_s (default 10)
%   'clust_method'    - clustering method (default 'kmeans')
%   'batch_method'    - 'harmony' | 'none' (default 'harmony')
%
%  --- Module 2 (sc_dock_ccc) ---
%   'lr_db'           - LR database table (default: lr_db_human.mat)
%   'condition'       - condition label per cell (default [])
%   'cond_labels'     - {cond1, cond2} for multi-group (default {})
%   'n_perm'          - permutation count (default 100)
%   'pval_cutoff'     - significance threshold (default 0.05)
%   'ccc_min_cells'   - min cells/type for CCC (default 10)
%   'sender'          - restrict sender cell types (default [])
%   'receiver'        - restrict receiver cell types (default [])
%
%  --- Module 3 (sc_dock_vina) ---
%   'compounds'       - compound IDs or struct array (default []; skip docking)
%   'proteins'        - protein IDs or struct array. If [], top targets from
%                       CCC result are extracted automatically. (default [])
%   'outdir'          - output directory (default './dock_output')
%   'vina_exe'        - path to Vina executable (default 'vina')
%   'obabel_exe'      - path to OpenBabel (default 'obabel')
%   'exhaustiveness'  - Vina exhaustiveness (default 8)
%   'n_modes'         - binding modes per run (default 9)
%   'global_docking'  - use whole-protein bounding box (default true)
%
% OUTPUT:
%   result  - struct with fields:
%     .rna    - output of sc_dock_rna
%     .ccc    - output of sc_dock_ccc  ([] if skipped)
%     .vina   - output of sc_dock_vina ([] if no compounds given)
%
% DEPENDENCIES:
%   sc_dock_rna, sc_dock_ccc, sc_dock_vina
%   scGEAToolbox_dev must be on the MATLAB path.
%
% See also: SC_DOCK_RNA, SC_DOCK_CCC, SC_DOCK_VINA

p = inputParser;
p.KeepUnmatched = true;
addRequired(p, 'X');
addRequired(p, 'g');
% Module 1 params
addParameter(p, 'batch_id',       [],         @(x) isempty(x) || isvector(x));
addParameter(p, 'species',        'human',    @ischar);
addParameter(p, 'libsz_cutoff',   500,        @isnumeric);
addParameter(p, 'mt_ratio',       0.15,       @isnumeric);
addParameter(p, 'min_cells',      10,         @isnumeric);
addParameter(p, 'n_hvg',          2000,       @isnumeric);
addParameter(p, 'n_pcs',          0,          @isnumeric);
addParameter(p, 'n_clusters',     10,         @isnumeric);
addParameter(p, 'clust_method',   'kmeans',   @ischar);
addParameter(p, 'batch_method',   'harmony',  @ischar);
% Module 2 params
addParameter(p, 'lr_db',          [],         @(x) isempty(x) || istable(x));
addParameter(p, 'condition',      [],         @(x) isempty(x) || isvector(x));
addParameter(p, 'cond_labels',    {},         @iscell);
addParameter(p, 'n_perm',         100,        @isnumeric);
addParameter(p, 'pval_cutoff',    0.05,       @isnumeric);
addParameter(p, 'ccc_min_cells',  10,         @isnumeric);
addParameter(p, 'sender',         [],         @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(p, 'receiver',       [],         @(x) isempty(x) || isstring(x) || iscellstr(x));
% Module 3 params
addParameter(p, 'compounds',      [],         @(x) true);
addParameter(p, 'proteins',       [],         @(x) true);
addParameter(p, 'outdir',         './dock_output', @ischar);
addParameter(p, 'vina_exe',       'vina',     @ischar);
addParameter(p, 'obabel_exe',     'obabel',   @ischar);
addParameter(p, 'exhaustiveness', 8,          @isnumeric);
addParameter(p, 'n_modes',        9,          @isnumeric);
addParameter(p, 'global_docking', true,       @islogical);
parse(p, X, g, varargin{:});
opt = p.Results;

result.rna  = [];
result.ccc  = [];
result.vina = [];

% =========================================================================
% MODULE 1: scRNA-seq preprocessing
% =========================================================================
fprintf('\n========== SC_DOCK: Module 1 — scRNA-seq Analysis ==========\n');
result.rna = sc_dock_rna(X, g, ...
    'batch_id',      opt.batch_id, ...
    'species',       opt.species, ...
    'libsz_cutoff',  opt.libsz_cutoff, ...
    'mt_ratio',      opt.mt_ratio, ...
    'min_cells',     opt.min_cells, ...
    'n_hvg',         opt.n_hvg, ...
    'n_pcs',         opt.n_pcs, ...
    'n_clusters',    opt.n_clusters, ...
    'clust_method',  opt.clust_method, ...
    'batch_method',  opt.batch_method);

% =========================================================================
% MODULE 2: Cell-cell communication
% =========================================================================
fprintf('\n========== SC_DOCK: Module 2 — Cell-Cell Communication ==========\n');
result.ccc = sc_dock_ccc(result.rna.X_norm, result.rna.g_hvg, ...
    result.rna.c_cell_type, ...
    'lr_db',       opt.lr_db, ...
    'condition',   opt.condition, ...
    'cond_labels', opt.cond_labels, ...
    'n_perm',      opt.n_perm, ...
    'pval_cutoff', opt.pval_cutoff, ...
    'min_cells',   opt.ccc_min_cells, ...
    'sender',      opt.sender, ...
    'receiver',    opt.receiver);

% =========================================================================
% MODULE 3: Molecular docking (only if compounds supplied)
% =========================================================================
if isempty(opt.compounds)
    fprintf('\n[sc_dock] No compounds provided — skipping docking module.\n');
    return
end

fprintf('\n========== SC_DOCK: Module 3 — Molecular Docking ==========\n');

% Auto-resolve protein targets from top CCC interactions if not provided
proteins = opt.proteins;
if isempty(proteins) && height(result.ccc.T_interactions) > 0
    proteins = sc_dock_gene2pdb(result.ccc.T_interactions, 'top_n', 20);
    fprintf('[sc_dock] Auto-selected %d target proteins from CCC results.\n', ...
        numel(proteins));
end

if isempty(proteins)
    warning('sc_dock:NoProteins', ...
        'No proteins specified and none could be inferred from CCC results. Skipping docking.');
    return
end

result.vina = sc_dock_vina(proteins, opt.compounds, ...
    'outdir',         opt.outdir, ...
    'vina_exe',       opt.vina_exe, ...
    'obabel_exe',     opt.obabel_exe, ...
    'exhaustiveness', opt.exhaustiveness, ...
    'n_modes',        opt.n_modes, ...
    'global_docking', opt.global_docking);

fprintf('\n========== SC_DOCK: Pipeline complete ==========\n');
fprintf('Top docking hits:\n');
if ~isempty(result.vina) && height(result.vina.T_results) > 0
    disp(result.vina.T_results(1:min(10, height(result.vina.T_results)), :));
end
end

% Gene->PDB mapping is handled by the public sc_dock_gene2pdb function.
