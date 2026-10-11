function results = run_sctenifoldxct(sample_id, celltype1, celltype2, data_dir, out_dir, null_control)
% LLM.RUN_SCTENIFOLDXCT  Ligand-receptor interaction analysis between two cell types.
%
%   results = llm.run_sctenifoldxct(sample_id, celltype1, celltype2)
%   results = llm.run_sctenifoldxct(sample_id, celltype1, celltype2, data_dir)
%   results = llm.run_sctenifoldxct(sample_id, celltype1, celltype2, data_dir, out_dir)
%   results = llm.run_sctenifoldxct(..., out_dir, null_control)
%
%   Detects ligand-receptor pairs co-expressed across two cell types within
%   a single GEO sample using spectral manifold alignment of gene regulatory
%   networks (scTenifoldXct). Runs both directions:
%     Direction 1: celltype1 (ligand source) -> celltype2 (receptor target)
%     Direction 2: celltype2 (ligand source) -> celltype1 (receptor target)
%
%   Method (Yang et al., Cell Systems 2023. PMID:36787742):
%     1. Build a PCNet GRN (principal-component regression) for each cell
%        type on log-normalized expression, via ten.i_pcnet
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
%     null_control - also run ten.xctnullcontrol for each direction
%                  (default true): 20 reruns with the networks' gene labels
%                  shuffled, about a minute in all at the runner's sizes.
%                  At the default settings scTenifoldXct's calls do not
%                  depend on the networks (pipelines/glyco_xct/Notes.md
%                  section 6k), and this is what tells a report so: each
%                  direction's verdict goes into the JSON as
%                  network_control, and the tables gain EmpP, ShuffleRatio,
%                  ShuffleSigFrac and ShuffleMedianDist. Pass false to skip.
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
%       .control1   - ten.xctnullcontrol's diagnostics for T1's direction
%                     (struct with Verdict, NetworkSensitive, MedianRho, ...),
%                     empty when null_control is false
%       .control2   - the same for T2's direction
%
%   Reference:
%     Yang Y et al. (2023). scTenifoldXct: A semi-supervised method for
%     predicting cell-cell interactions and mapping cellular communication
%     graphs. Cell Systems 14(4):302-311.e4. PMID:36787742.
%
%   Example (from agent via evaluate_matlab_code):
%     results = llm.run_sctenifoldxct('GSM2333580', 'T cell', 'Cancer cell', ...
%                   'C:/abs/path/to/data', 'C:/abs/path/to/output');

if nargin < 4 || isempty(data_dir), data_dir = 'data'; end
if nargin < 5, out_dir = []; end
if nargin < 6 || isempty(null_control), null_control = true; end

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
    'T2sig', table(), ...
    'control1', struct([]), ...
    'control2', struct([]));

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
sce = sce.qcfilterwhitelist(500, 0.15, 15, 200, []);
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
fprintf('%s: %s\n', celltype1, pkg.i_plural(n1, 'cell'));
fprintf('%s: %s\n', celltype2, pkg.i_plural(n2, 'cell'));

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

T_pair = [];
try
    % Returns {T1, T2}: T1 = ct1->ct2, T2 = ct2->ct1
    [T_pair, grns] = ten.sctenifoldxct(sce, celltype1, celltype2, true);
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

% ---- Shuffled-network control ---------------------------------------
% Each direction's control re-runs it with the networks the run above
% built (swapped for the reverse direction) and returns the same table with
% the control's columns, so it replaces T1/T2. Seed 1 keeps reports
% reproducible. A failure here costs the control, not the result.
if null_control
    nets = {grns, struct('A_s', grns.A_t, 'A_t', grns.A_s, 'genes', grns.genes)};
    pairs = [celltype1, celltype2; celltype2, celltype1];
    Ts = {T1, T2};
    controls = {struct([]), struct([])};
    for d = 1:2
        if height(Ts{d}) < 3
            continue   % too few candidates to compare rankings
        end
        fprintf('Shuffled-network control, %s -> %s ...\n', pairs(d, 1), pairs(d, 2));
        try
            [Tn, Sn] = ten.xctnullcontrol(sce, pairs(d, 1), pairs(d, 2), ...
                Grns=nets{d}, Seed=1, Verbose=false);
            fprintf('  %s\n', Sn.Verdict);
            Ts{d} = sortrows(Tn, 'dist');
            controls{d} = Sn;
        catch ME
            warning('llm:run_sctenifoldxct:controlFailed', ...
                'Shuffled-network control failed for %s -> %s: %s', ...
                pairs(d, 1), pairs(d, 2), ME.message);
        end
    end
    [T1, T2] = Ts{:};
    [results.control1, results.control2] = controls{:};
end

% Significant pairs (p_value < 0.05)
T1sig = i_significant(T1);
T2sig = i_significant(T2);

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


% ---- Helper: pairs with p_value < 0.05 ------------------------------
function Tsig = i_significant(T)
% A table without a p_value column has nothing to call significant. This
% used to return the whole table, so every candidate pair went to the
% report -- and to the agent's interpretation -- as a significant one.
if isempty(T)
    Tsig = T;
elseif ismember('p_value', T.Properties.VariableNames)
    Tsig = T(T.p_value < 0.05, :);
else
    warning('llm:run_sctenifoldxct:noPValue', ...
        'The scTenifoldXct result has no p_value column; no pair is reported as significant.');
    Tsig = T([], :);
end
end


% ---- Helper: print JSON summary -------------------------------------
function i_print_json_summary(results, top_n)
if nargin < 2, top_n = 20; end

dir1 = struct( ...
    'ligand_cell',    results.celltype1, ...
    'receptor_cell',  results.celltype2, ...
    'n_total_pairs',  height(results.T1), ...
    'n_sig_pairs',    height(results.T1sig), ...
    'top_pairs',      {i_table_to_structs(results.T1sig, top_n)}, ...
    'network_control', i_control_json(results.control1));

dir2 = struct( ...
    'ligand_cell',    results.celltype2, ...
    'receptor_cell',  results.celltype1, ...
    'n_total_pairs',  height(results.T2), ...
    'n_sig_pairs',    height(results.T2sig), ...
    'top_pairs',      {i_table_to_structs(results.T2sig, top_n)}, ...
    'network_control', i_control_json(results.control2));

% STATUS for the schema GEOcellar's analysis node reads. A run that gets
% here has results; zero significant pairs is one of them, reported by
% n_sig_pairs, and an analysis that cannot run raises an error instead.
summary = struct( ...
    'status',             'completed', ...
    'reason',             '', ...
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
hasControl = ismember('EmpP', T.Properties.VariableNames);
for i = 1:n
    s = struct('ligand',   char(T.ligand(i)), ...
               'receptor', char(T.receptor(i)), ...
               'dist',     round(double(T.dist(i)), 6), ...
               'p_value',  double(T.p_value(i)));
    if hasControl
        % The pair against its own shuffled networks; see i_control_json
        s.emp_p = round(double(T.EmpP(i)), 4);
        s.shuffle_ratio = round(double(T.ShuffleRatio(i)), 4);
    end
    rows{end+1} = s; %#ok<AGROW>
end
end

% ---- Helper: shuffled-network control for JSON ----------------------
function c = i_control_json(S)
% ten.xctnullcontrol's verdict for one direction, or a note that it did
% not run. VERDICT is a full sentence meant to be read as it stands.
if isempty(S)
    c = struct('ran', false, 'verdict', '', 'network_sensitive', [], ...
        'median_spearman_vs_shuffled', [], 'median_sig_jaccard_vs_shuffled', [], ...
        'n_pairs_beating_shuffles', [], 'n_shuffles', 0);
    return
end
c = struct('ran', true, 'verdict', char(S.Verdict), ...
    'network_sensitive', S.NetworkSensitive, ...
    'median_spearman_vs_shuffled', round(S.MedianRho, 4), ...
    'median_sig_jaccard_vs_shuffled', round(S.MedianSigJaccard, 4), ...
    'n_pairs_beating_shuffles', S.NumEmpPBelow05, ...
    'n_shuffles', S.NumShuffles);
end


% ---- Helper: locate and load cleandata.mat --------------------------
function sce = i_load_sce(sample_id, data_dir)
% The shared loader. This file carried its own copy of it.
sce = llm.i_load_sce(sample_id, data_dir);
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
% A local seeded stream: identical calls pick the same cells, and the
% caller's global stream is left alone.
stream = RandStream("twister", "Seed", 1);
for ctname = [ct1, ct2]
    mask = ct == ctname;
    idx  = find(mask);
    if numel(idx) > max_cells
        keep = idx(randperm(stream, numel(idx), max_cells));
        drop = setdiff(idx, keep);
        keep_all = setdiff(1:sce.NumCells, drop);
        sce = sce.selectcells(keep_all);
        ct  = sce.c_cell_type_tx;
        fprintf('  Subsampled "%s": %d → %d cells\n', ctname, numel(idx), max_cells);
    end
end
end


