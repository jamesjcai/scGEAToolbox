function results = run_sctenifoldxct(sample_id, celltype1, celltype2, data_dir, out_dir)
% LLM.RUN_SCTENIFOLDXCT  Ligand-receptor interaction analysis between two cell types.
%
%   results = llm.run_sctenifoldxct(sample_id, celltype1, celltype2)
%   results = llm.run_sctenifoldxct(sample_id, celltype1, celltype2, data_dir)
%   results = llm.run_sctenifoldxct(sample_id, celltype1, celltype2, data_dir, out_dir)
%
%   Detects ligand-receptor pairs co-expressed across two cell types within
%   a single GEO sample using spectral manifold alignment of gene regulatory
%   networks (scTenifoldXct). Runs both directions:
%     Direction 1: celltype1 (ligand source) -> celltype2 (receptor target)
%     Direction 2: celltype2 (ligand source) -> celltype1 (receptor target)
%
%   Method (Ma et al., Cell Systems 2023. PMID:36787742):
%     1. Build partial-correlation GRNs for each cell type via net.pcrnet
%     2. Assemble block weight matrix with L-R database correspondences
%     3. Spectral manifold alignment via graph Laplacian eigenvectors
%     4. Rank L-R pairs by embedding distance (small dist = strong interaction)
%     5. Nonparametric left-tail null test against random gene pairs
%
%   Inputs:
%     sample_id  - GSM accession string (e.g. 'GSM2333580')
%     celltype1  - first cell type label (string, must match sce.c_cell_type_tx)
%     celltype2  - second cell type label (string, must match sce.c_cell_type_tx)
%     data_dir   - root folder containing downloaded .mat files, organised as
%                  data_dir/<study_id>/<sample_id>/cleandata.mat
%                  Default: 'data' (relative to MATLAB working directory).
%                  Pass an absolute path when calling from the agent.
%     out_dir    - folder to write Excel result files (optional).
%                  If omitted or empty, no files are written.
%
%   Output:
%     results - struct with fields:
%       .sample     - GSM accession
%       .celltype1  - first cell type label
%       .celltype2  - second cell type label
%       .n1         - number of cells in celltype1
%       .n2         - number of cells in celltype2
%       .T1         - L-R table, celltype1->celltype2 (all pairs, sorted by dist)
%                     Columns: ligand, receptor, dist, correspondence, p_value
%       .T2         - L-R table, celltype2->celltype1 (all pairs, sorted by dist)
%       .T1sig      - significant pairs from T1 (p_value < 0.05)
%       .T2sig      - significant pairs from T2 (p_value < 0.05)
%
%   Reference:
%     Ma W et al. (2023). Characterizing cell-cell communication in
%     single-cell multi-omics data with scTenifoldXct.
%     Cell Systems 16(4):337-354. PMID:36787742.
%
%   Example (from agent via evaluate_matlab_code):
%     results = llm.run_sctenifoldxct('GSM2333580', 'T cell', 'Cancer cell', ...
%                   'C:/abs/path/to/data', 'C:/abs/path/to/output');

if nargin < 4 || isempty(data_dir), data_dir = 'data'; end
if nargin < 5, out_dir = []; end

celltype1 = string(celltype1);
celltype2 = string(celltype2);

results = struct( ...
    'sample',    char(sample_id), ...
    'celltype1', char(celltype1), ...
    'celltype2', char(celltype2), ...
    'n1',    0, ...
    'n2',    0, ...
    'T1',    table(), ...
    'T2',    table(), ...
    'T1sig', table(), ...
    'T2sig', table());

% ---- Load SCE object ------------------------------------------------
sce = i_load_sce(sample_id, data_dir);
fprintf('Sample (%s): %d genes x %d cells\n', sample_id, sce.NumGenes, sce.NumCells);

% ---- Limit to top 3000 HVGs for speed -------------------------------
max_genes = 3000;
if sce.NumGenes > max_genes
    fprintf('Limiting to top %d HVGs (from %d genes)...\n', max_genes, sce.NumGenes);
    sce = llm.i_limit_genes(sce, max_genes);
    fprintf('Genes after HVG selection: %d\n', sce.NumGenes);
end

% ---- QC filter ------------------------------------------------------
fprintf('Applying QC filter...\n');
sce = sce.qcfilterwhitelist(1000, 0.15, 15, 500, []);
fprintf('After QC: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);

% ---- List available cell types --------------------------------------
ct = sce.c_cell_type_tx;
available_ct = unique(ct);
available_ct = available_ct(~strcmpi(available_ct, 'undetermined'));
fprintf('Available cell types (%d): %s\n', numel(available_ct), strjoin(available_ct, ', '));

% ---- Validate requested cell types ----------------------------------
if ~any(available_ct == celltype1)
    error('llm:run_sctenifoldxct:celltypeNotFound', ...
        'Cell type "%s" not found in sample %s.\nAvailable: %s', ...
        celltype1, sample_id, strjoin(available_ct, ', '));
end
if ~any(available_ct == celltype2)
    error('llm:run_sctenifoldxct:celltypeNotFound', ...
        'Cell type "%s" not found in sample %s.\nAvailable: %s', ...
        celltype2, sample_id, strjoin(available_ct, ', '));
end

max_cells = 1000;   % subsample each cell type; GRN build time grows with cells

n1 = sum(ct == celltype1);
n2 = sum(ct == celltype2);
results.n1 = n1;
results.n2 = n2;
fprintf('%s: %d cells\n', celltype1, n1);
fprintf('%s: %d cells\n', celltype2, n2);

% Subsample cell types that exceed the cap
sce = i_subsample_celltypes(sce, ct, celltype1, celltype2, max_cells);

if n1 < 50
    warning('llm:run_sctenifoldxct:fewCells', ...
        '"%s" has only %d cells; GRN estimation may be unreliable (<50).', celltype1, n1);
end
if n2 < 50
    warning('llm:run_sctenifoldxct:fewCells', ...
        '"%s" has only %d cells; GRN estimation may be unreliable (<50).', celltype2, n2);
end

% ---- Prepare output directory ---------------------------------------
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

i_log(out_dir, sprintf('START scTenifoldXct: %s | %s (n=%d) <-> %s (n=%d)', ...
    sample_id, celltype1, n1, celltype2, n2));

% ---- Run scTenifoldXct (both directions) ----------------------------
fprintf('\nRunning scTenifoldXct: %s <-> %s ...\n', celltype1, celltype2);
fprintf('WARNING: GRN construction may take several minutes.\n');

T_pair = [];
try
    % Returns {T1, T2}: T1 = ct1->ct2, T2 = ct2->ct1
    T_pair = ten.sctenifoldxct(sce, celltype1, celltype2, true);
catch ME
    error('llm:run_sctenifoldxct:failed', ...
        'scTenifoldXct failed: %s', ME.message);
end

% Unpack results
if iscell(T_pair) && numel(T_pair) == 2
    T1 = T_pair{1};
    T2 = T_pair{2};
else
    T1 = T_pair;
    T2 = table();
end

if ~istable(T1), T1 = table(); end
if ~istable(T2), T2 = table(); end

% Significant pairs (p_value < 0.05)
if ~isempty(T1) && ismember('p_value', T1.Properties.VariableNames)
    T1sig = T1(T1.p_value < 0.05, :);
else
    T1sig = T1;
end
if ~isempty(T2) && ismember('p_value', T2.Properties.VariableNames)
    T2sig = T2(T2.p_value < 0.05, :);
else
    T2sig = T2;
end

fprintf('%s -> %s: %d total pairs, %d significant (p<0.05)\n', ...
    celltype1, celltype2, height(T1), height(T1sig));
fprintf('%s -> %s: %d total pairs, %d significant (p<0.05)\n', ...
    celltype2, celltype1, height(T2), height(T2sig));

i_log(out_dir, sprintf('DONE: %s->%s %d sig pairs; %s->%s %d sig pairs', ...
    celltype1, celltype2, height(T1sig), celltype2, celltype1, height(T2sig)));

results.T1    = T1;
results.T2    = T2;
results.T1sig = T1sig;
results.T2sig = T2sig;

% ---- Save Excel files -----------------------------------------------
if ~isempty(out_dir)
    ct1_safe = matlab.lang.makeValidName(char(celltype1));
    ct2_safe = matlab.lang.makeValidName(char(celltype2));
    gsm_safe = matlab.lang.makeValidName(char(sample_id));

    % Direction 1
    if ~isempty(T1sig)
        f1 = fullfile(out_dir, sprintf('XCT_%s_%s_to_%s.xlsx', gsm_safe, ct1_safe, ct2_safe));
        try
            writetable(T1sig, f1, 'FileType', 'spreadsheet', 'Sheet', 'Significant_pairs');
            writetable(T1,    f1, 'FileType', 'spreadsheet', 'Sheet', 'All_pairs');
            fprintf('Saved: %s\n', f1);
        catch ME
            warning('Could not save %s: %s', f1, ME.message);
        end
    end

    % Direction 2
    if ~isempty(T2sig)
        f2 = fullfile(out_dir, sprintf('XCT_%s_%s_to_%s.xlsx', gsm_safe, ct2_safe, ct1_safe));
        try
            writetable(T2sig, f2, 'FileType', 'spreadsheet', 'Sheet', 'Significant_pairs');
            writetable(T2,    f2, 'FileType', 'spreadsheet', 'Sheet', 'All_pairs');
            fprintf('Saved: %s\n', f2);
        catch ME
            warning('Could not save %s: %s', f2, ME.message);
        end
    end
end

fprintf('\nscTenifoldXct complete.\n');
i_log(out_dir, 'COMPLETE scTenifoldXct: emitting JSON summary');

% ---- Print JSON summary for agent consumption -----------------------
i_print_json_summary(results);
end


% ---- Helper: print JSON summary -------------------------------------
function i_print_json_summary(results, top_n)
if nargin < 2, top_n = 20; end

dir1 = struct( ...
    'ligand_cell',    results.celltype1, ...
    'receptor_cell',  results.celltype2, ...
    'n_total_pairs',  height(results.T1), ...
    'n_sig_pairs',    height(results.T1sig), ...
    'top_pairs',      {i_table_to_structs(results.T1sig, top_n)});

dir2 = struct( ...
    'ligand_cell',    results.celltype2, ...
    'receptor_cell',  results.celltype1, ...
    'n_total_pairs',  height(results.T2), ...
    'n_sig_pairs',    height(results.T2sig), ...
    'top_pairs',      {i_table_to_structs(results.T2sig, top_n)});

summary = struct( ...
    'sample',             results.sample, ...
    'celltype1',          results.celltype1, ...
    'celltype2',          results.celltype2, ...
    'n_cells_ct1',        results.n1, ...
    'n_cells_ct2',        results.n2, ...
    'direction_ct1_to_ct2', dir1, ...
    'direction_ct2_to_ct1', dir2);

fprintf('\n%%SCTENIFOLDXCT_JSON_SUMMARY_BEGIN%%\n%s\n%%SCTENIFOLDXCT_JSON_SUMMARY_END%%\n', ...
    jsonencode(summary, 'PrettyPrint', true));
end


% ---- Helper: extract top rows as struct list for JSON ---------------
function rows = i_table_to_structs(T, n)
rows = {};
if isempty(T) || ~istable(T), return; end
n = min(n, height(T));
for i = 1:n
    s = struct('ligand',   char(T.ligand(i)), ...
               'receptor', char(T.receptor(i)), ...
               'dist',     round(double(T.dist(i)), 6), ...
               'p_value',  double(T.p_value(i)));
    rows{end+1} = s; %#ok<AGROW>
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
        error('llm:run_sctenifoldxct:fileNotFound', ...
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


% ---- Helper: append a timestamped line to progress.log --------------
function i_log(out_dir, msg)
if isempty(out_dir), return; end
try
    fid = fopen(fullfile(out_dir, 'progress.log'), 'a');
    if fid ~= -1
        fprintf(fid, '[%s] %s\n', datestr(now, 'yyyy-mm-dd HH:MM:SS'), msg); %#ok<TNOW1,DATST>
        fclose(fid);
    end
catch
    % progress log is best-effort; never fail the main task
end
end


% ---- Helper: subsample two named cell types to at most max_cells each
function sce = i_subsample_celltypes(sce, ct, ct1, ct2, max_cells)
for ctname = [ct1, ct2]
    mask = ct == ctname;
    idx  = find(mask);
    if numel(idx) > max_cells
        keep = idx(randperm(numel(idx), max_cells));
        drop = setdiff(idx, keep);
        keep_all = setdiff(1:sce.NumCells, drop);
        sce = sce.selectcells(keep_all);
        ct  = sce.c_cell_type_tx;
        fprintf('  Subsampled "%s": %d → %d cells\n', ctname, numel(idx), max_cells);
    end
end
end


