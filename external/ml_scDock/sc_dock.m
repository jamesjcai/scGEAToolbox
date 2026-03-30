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
    proteins = i_extract_top_proteins(result.ccc.T_interactions, opt.species);
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

% =========================================================================
% Local helper: infer protein PDB IDs from top-ranked CCC interactions
% =========================================================================
function proteins = i_extract_top_proteins(T_int, species) %#ok<INUSD>
% Collect unique receptor gene names from top 20 interactions,
% then look up PDB IDs from the bundled receptor_reference.csv.
% Falls back to UniProt web query for genes not in the reference.
top_n     = min(20, height(T_int));
rec_genes = unique(T_int.receptor(1:top_n));

% Load receptor reference (gene -> PDB chain, from scDock repo)
ref_map = i_load_receptor_ref();

proteins = strings(0);
for i = 1:numel(rec_genes)
    gene = char(rec_genes(i));
    if isKey(ref_map, gene)
        % May contain multiple models separated by "; " — take the first
        entries = strsplit(ref_map(gene), '; ');
        pdb_chain = strtrim(entries{1});   % e.g. "1A22.B"
        pdb_id    = extractBefore(pdb_chain, '.');  % e.g. "1A22"
        if ~isempty(pdb_id)
            proteins(end+1) = string(pdb_id); %#ok<AGROW>
        end
    else
        pdb_id = i_gene2pdb_uniprot(gene);
        if ~isempty(pdb_id)
            proteins(end+1) = string(pdb_id); %#ok<AGROW>
        end
    end
end
proteins = unique(proteins);
end

function ref_map = i_load_receptor_ref()
% Load receptor_reference.csv (protein_name, PDB_model) into a Map.
% Searches same directory as this file, then current directory.
ref_map = containers.Map('KeyType','char','ValueType','char');
candidates = {
    fullfile(fileparts(mfilename('fullpath')), 'receptor_reference.csv');
    fullfile(pwd, 'receptor_reference.csv');
};
f = '';
for k = 1:numel(candidates)
    if isfile(candidates{k}), f = candidates{k}; break; end
end
if isempty(f)
    warning('sc_dock:NoReceptorRef', ...
        'receptor_reference.csv not found. Falling back to UniProt lookup.');
    return
end
T = readtable(f, 'TextType', 'string');
for k = 1:height(T)
    ref_map(char(T.protein_name(k))) = char(T.PDB_model(k));
end
end

function pdb_id = i_gene2pdb_uniprot(gene)
% Fallback: query UniProt REST API for gene -> best PDB structure.
pdb_id = '';
try
    url = sprintf( ...
        'https://rest.uniprot.org/uniprotkb/search?query=gene:%s+organism_id:9606+reviewed:true&format=tsv&fields=accession,xref_pdb', ...
        urlencode(gene));
    txt  = webread(url);
    lines = strsplit(strtrim(txt), newline);
    if numel(lines) < 2, return; end
    cols = strsplit(lines{2}, char(9));
    if numel(cols) >= 2 && ~isempty(strtrim(cols{2}))
        pdbs   = strsplit(strtrim(cols{2}), ';');
        pdb_id = strtrim(pdbs{1});
    end
catch
end
end
