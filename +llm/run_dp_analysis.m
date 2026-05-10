function results = run_dp_analysis(sample_id1, sample_id2, data_dir, out_dir, gene_set_option, species)
% LLM.RUN_DP_ANALYSIS  Cell type-specific differential program analysis between two GEO samples.
%
%   results = llm.run_dp_analysis(sample_id1, sample_id2)
%   results = llm.run_dp_analysis(sample_id1, sample_id2, data_dir)
%   results = llm.run_dp_analysis(sample_id1, sample_id2, data_dir, out_dir)
%   results = llm.run_dp_analysis(sample_id1, sample_id2, data_dir, out_dir, gene_set_option)
%   results = llm.run_dp_analysis(sample_id1, sample_id2, data_dir, out_dir, gene_set_option, species)
%
%   Loads cleandata.mat for each sample (each file contains a
%   SingleCellExperiment variable named 'sce') and performs differential
%   program (DP) analysis using sc_dpg for every cell type shared between
%   the two samples.
%
%   The non-interactive gene-set collections supported here are:
%     2 / 'TF'          - DoRothEA TF target programs
%     3 / 'Predefined'  - scGEAToolbox predefined marker programs
%
%   MSigDB option 1 is intentionally rejected here because pkg.e_getgenesets(1)
%   opens an interactive selector dialog, which is not suitable for agent use.
%
%   If out_dir is provided, results are saved as Excel files
%   (DP_<id1>_vs_<id2>_<celltype>.xlsx) with sheets: All programs,
%   Up-regulated, Down-regulated.
%
%   Inputs:
%     sample_id1       - GSM accession of sample 1 (e.g. 'GSM2333580')
%     sample_id2       - GSM accession of sample 2 (e.g. 'GSM2333581')
%     data_dir         - root folder containing downloaded .mat files, organised
%                        as data_dir/<study_id>/<sample_id>/cleandata.mat
%                        Default: 'data' (relative to MATLAB working directory).
%                        Pass an absolute path when calling from the agent.
%     out_dir          - folder to write Excel result files (optional).
%                        If omitted or empty, no files are written.
%     gene_set_option  - gene set collection for pkg.e_getgenesets.
%                        Default: 2 ('TF').
%     species          - species passed to pkg.e_getgenesets. Default: 'human'.
%
%   Output:
%     results - struct array with one element per shared cell type:
%       .cell_type       - cell type label (string)
%       .n1              - number of cells from sample 1
%       .n2              - number of cells from sample 2
%       .gene_set_option - collection identifier used for the run
%       .T               - full DP table
%       .Tup             - programs higher in sample 1
%       .Tdn             - programs higher in sample 2
%
%   Example (from agent via evaluate_matlab_code):
%     results = llm.run_dp_analysis('GSM2333580', 'GSM2333581', ...
%                   'C:/abs/path/to/data', 'C:/abs/path/to/output', 2, 'human');

if nargin < 3 || isempty(data_dir), data_dir = 'data'; end
if nargin < 4, out_dir = []; end
if nargin < 5 || isempty(gene_set_option), gene_set_option = 2; end
if nargin < 6 || isempty(species), species = 'human'; end

max_cells = 2000;   % subsample per cell type for speed

i_validate_gene_set_option(gene_set_option);

results = struct('cell_type', {}, 'n1', {}, 'n2', {}, ...
    'gene_set_option', {}, 'T', {}, 'Tup', {}, 'Tdn', {});

% ---- Load SCE objects -----------------------------------------------
sce1 = i_load_sce(sample_id1, data_dir);
sce2 = i_load_sce(sample_id2, data_dir);

fprintf('Sample 1 (%s): %d genes x %d cells\n', sample_id1, sce1.NumGenes, sce1.NumCells);
fprintf('Sample 2 (%s): %d genes x %d cells\n', sample_id2, sce2.NumGenes, sce2.NumCells);

% ---- Load gene sets --------------------------------------------------
[setmatrx, setnames, setgenes] = pkg.e_getgenesets(gene_set_option, species);
if isempty(setmatrx) || isempty(setnames) || isempty(setgenes)
    error('llm:run_dp_analysis:noGeneSets', ...
        'No gene sets were loaded for option "%s".', string(gene_set_option));
end
fprintf('Gene set collection "%s": %d programs across %d genes\n', ...
    string(gene_set_option), numel(setnames), numel(setgenes));

% ---- Align to common gene set ---------------------------------------
[common_genes, idx1, idx2] = intersect(sce1.g, sce2.g, 'stable');
if isempty(common_genes)
    error('llm:run_dp_analysis:noCommonGenes', ...
        'No common genes between %s and %s.', sample_id1, sample_id2);
end
fprintf('Common genes: %d\n', numel(common_genes));

X1_all = log1p(sc_norm(sce1.X(idx1, :)));
X2_all = log1p(sc_norm(sce2.X(idx2, :)));

% ---- Identify shared cell types -------------------------------------
ct1 = sce1.c_cell_type_tx;
ct2 = sce2.c_cell_type_tx;

shared_ct = intersect(unique(ct1), unique(ct2));
shared_ct = shared_ct(~strcmpi(shared_ct, 'undetermined'));

if isempty(shared_ct)
    warning('llm:run_dp_analysis:noCellTypeAnnotation', ...
        'No shared annotated cell types found. Running DP on all cells combined.');
    shared_ct = "all_cells";
    ct1(:) = "all_cells";
    ct2(:) = "all_cells";
end

fprintf('Shared cell types (%d): %s\n', numel(shared_ct), strjoin(shared_ct, ', '));

% ---- Prepare output directory ---------------------------------------
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

% ---- DP per cell type -----------------------------------------------
ri = 0;
for k = 1:numel(shared_ct)
    ct = shared_ct(k);
    mask1 = ct1 == ct;
    mask2 = ct2 == ct;
    n1 = sum(mask1);
    n2 = sum(mask2);

    if n1 < 500 || n2 < 500
        fprintf('Skipping "%s": fewer than 500 cells (%d in sample1, %d in sample2).\n', ...
            ct, n1, n2);
        continue;
    end

    fprintf('DP for "%s": %d vs %d cells ... ', ct, n1, n2);

    idx1 = find(mask1);
    idx2 = find(mask2);
    if numel(idx1) > max_cells
        idx1 = idx1(randperm(numel(idx1), max_cells));
        fprintf('(subsampled→%d) ', max_cells);
    end
    if numel(idx2) > max_cells
        idx2 = idx2(randperm(numel(idx2), max_cells));
        fprintf('(subsampled→%d) ', max_cells);
    end
    n1 = numel(idx1);
    n2 = numel(idx2);

    try
        T = sc_dpg(X1_all(:, idx1), X2_all(:, idx2), common_genes, ...
            setmatrx, setnames, setgenes);
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        continue;
    end

    Tup = T(T.avg_log2FC > 0, :);
    Tdn = T(T.avg_log2FC < 0, :);

    fprintf('done. Up: %d  Down: %d  (FDR<0.01, |log2FC|>=1)\n', ...
        height(Tup), height(Tdn));

    ri = ri + 1;
    results(ri).cell_type       = ct;
    results(ri).n1              = n1;
    results(ri).n2              = n2;
    results(ri).gene_set_option = string(gene_set_option);
    results(ri).T               = T;
    results(ri).Tup             = Tup;
    results(ri).Tdn             = Tdn;

    if ~isempty(out_dir)
        outfile = sprintf('DP_%s_vs_%s_%s.xlsx', ...
            matlab.lang.makeValidName(sample_id1), ...
            matlab.lang.makeValidName(sample_id2), ...
            matlab.lang.makeValidName(string(ct)));
        filesaved = fullfile(out_dir, outfile);
        try
            writetable(T, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All programs');
            writetable(Tup, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Down-regulated');
            fprintf('  Saved: %s\n', filesaved);
        catch ME
            warning('Could not save %s: %s', filesaved, ME.message);
        end
    end
end

fprintf('\nDP analysis complete: %d cell type(s) analysed.\n', numel(results));

% ---- Print JSON summary for agent consumption -----------------------
i_print_json_summary(results, sample_id1, sample_id2, gene_set_option);
end


function i_validate_gene_set_option(gene_set_option)
if isequal(gene_set_option, 1) || strcmpi(string(gene_set_option), "MSIGDB") || ...
        strcmpi(string(gene_set_option), "MSigDB Molecular Signatures")
    error('llm:run_dp_analysis:interactiveGeneSetOption', ...
        ['MSigDB gene set selection is interactive and not supported by ' ...
        'llm.run_dp_analysis. Use gene_set_option 2 (TF) or 3 (Predefined).']);
end
end


% ---- Helper: print JSON summary of top differential programs --------
function i_print_json_summary(results, sample_id1, sample_id2, gene_set_option, top_n)
if nargin < 5, top_n = 20; end

cell_types = {};
for k = 1:numel(results)
    r = results(k);

    top_up = i_table_to_structs(r.Tup, top_n);
    top_dn = i_table_to_structs(r.Tdn, top_n);

    cell_types{end+1} = struct( ...
        'cell_type', char(r.cell_type), ...
        'n1',        r.n1, ...
        'n2',        r.n2, ...
        'n_up',      height(r.Tup), ...
        'n_dn',      height(r.Tdn), ...
        'top_up',    {top_up}, ...
        'top_dn',    {top_dn}); %#ok<AGROW>
end

summary = struct( ...
    'sample1',          char(sample_id1), ...
    'sample2',          char(sample_id2), ...
    'gene_set_option',  char(string(gene_set_option)), ...
    'cell_types',       {cell_types});

fprintf('\n%%DP_JSON_SUMMARY_BEGIN%%\n%s\n%%DP_JSON_SUMMARY_END%%\n', ...
    jsonencode(summary, 'PrettyPrint', true));
end


% ---- Helper: extract top rows as struct list for JSON ---------------
function rows = i_table_to_structs(T, n)
rows = {};
if isempty(T), return; end
n = min(n, height(T));
for i = 1:n
    rows{end+1} = struct( ...
        'program',     char(T.setnames(i)), ...
        'log2FC',      round(T.avg_log2FC(i), 3), ...
        'p_val_adj',   T.p_val_adj(i), ...
        'gsetsize',    T.gsetsize(i)); %#ok<AGROW>
end
end


% ---- Helper: locate and load cleandata.mat --------------------------
function sce = i_load_sce(sample_id, data_dir)
hits = dir(fullfile(data_dir, '*', sample_id, 'cleandata.mat'));
if isempty(hits)
    flat = fullfile(data_dir, sample_id, 'cleandata.mat');
    if isfile(flat)
        mat_path = flat;
    else
        error('llm:run_dp_analysis:fileNotFound', ...
            'Cannot find cleandata.mat for sample "%s" under "%s".', ...
            sample_id, data_dir);
    end
else
    mat_path = fullfile(hits(1).folder, hits(1).name);
end
fprintf('Loading %s\n', mat_path);
s = load(mat_path, 'sce');
sce = s.sce;
end
