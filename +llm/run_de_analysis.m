function results = run_de_analysis(sample_id1, sample_id2, data_dir, out_dir)
% LLM.RUN_DE_ANALYSIS  Cell type-specific DE analysis between two GEO samples.
%
%   results = llm.run_de_analysis(sample_id1, sample_id2)
%   results = llm.run_de_analysis(sample_id1, sample_id2, data_dir)
%   results = llm.run_de_analysis(sample_id1, sample_id2, data_dir, out_dir)
%
%   Loads cleandata.mat for each sample (each file contains a
%   SingleCellExperiment variable named 'sce') and performs DE analysis
%   (Mann-Whitney U test via sc_deg) for every cell type shared between
%   the two samples. Post-processing follows the same pipeline as the
%   interactive batch DE tool (callback_DEGene2GroupsBatch):
%     - pkg.in_DETableProcess  : label columns with sample IDs, fix Inf
%     - pkg.e_processdetable   : split into Tup/Tdn with standard thresholds
%                                (min Δpct 5%, min |log2FC| 1.0, adj-p < 0.01)
%   If out_dir is provided, results are saved as Excel files
%   (DE_<id1>_vs_<id2>_<celltype>.xlsx) with sheets: All_genes,
%   Up-regulated, Down-regulated, Note.
%
%   Inputs:
%     sample_id1 - GSM accession of sample 1 (e.g. 'GSM2333580')
%     sample_id2 - GSM accession of sample 2 (e.g. 'GSM2333581')
%     data_dir   - root folder containing downloaded .mat files, organised
%                  as data_dir/<study_id>/<sample_id>/cleandata.mat
%                  Default: 'data' (relative to MATLAB working directory).
%                  Pass an absolute path when calling from the agent.
%     out_dir    - folder to write Excel result files (optional).
%                  If omitted or empty, no files are written.
%
%   Output:
%     results - struct array with one element per shared cell type:
%       .cell_type  - cell type label (string)
%       .n1         - number of cells from sample 1
%       .n2         - number of cells from sample 2
%       .T          - full DE table (all genes, labelled columns)
%       .Tup        - up-regulated DEGs (higher in sample 1)
%       .Tdn        - down-regulated DEGs (higher in sample 2)
%
%   Example (from agent via evaluate_matlab_code):
%     results = llm.run_de_analysis('GSM2333580', 'GSM2333581', ...
%                   'C:/abs/path/to/data', 'C:/abs/path/to/output');
%     for k = 1:numel(results)
%         fprintf('=== %s ===  up: %d  dn: %d\n', ...
%             results(k).cell_type, height(results(k).Tup), height(results(k).Tdn));
%         disp(results(k).Tup(1:min(10,end), :));
%     end

if nargin < 3 || isempty(data_dir), data_dir = 'data'; end
if nargin < 4, out_dir = []; end

% Default filter parameters matching gui.i_degparamset(true) defaults:
%   mindiffpct=0.05, minabsolfc=1.0, apvaluecut=0.01, sort by adj p-value
paramset = {0.05, 1.0, 0.01, 'Adjusted P-value'};

results = struct('cell_type', {}, 'n1', {}, 'n2', {}, ...
                 'T', {}, 'Tup', {}, 'Tdn', {});

% ---- Load SCE objects -----------------------------------------------
sce1 = i_load_sce(sample_id1, data_dir);
sce2 = i_load_sce(sample_id2, data_dir);

fprintf('Sample 1 (%s): %d genes x %d cells\n', sample_id1, sce1.NumGenes, sce1.NumCells);
fprintf('Sample 2 (%s): %d genes x %d cells\n', sample_id2, sce2.NumGenes, sce2.NumCells);

% ---- Align to common gene set ---------------------------------------
[common_genes, idx1, idx2] = intersect(sce1.g, sce2.g, 'stable');
if isempty(common_genes)
    error('llm:run_de_analysis:noCommonGenes', ...
        'No common genes between %s and %s.', sample_id1, sample_id2);
end
fprintf('Common genes: %d\n', numel(common_genes));

X1_all = sce1.X(idx1, :);
X2_all = sce2.X(idx2, :);

% ---- Identify shared cell types -------------------------------------
ct1 = sce1.c_cell_type_tx;
ct2 = sce2.c_cell_type_tx;

shared_ct = intersect(unique(ct1), unique(ct2));
shared_ct = shared_ct(~strcmpi(shared_ct, 'undetermined'));

if isempty(shared_ct)
    warning('llm:run_de_analysis:noCellTypeAnnotation', ...
        'No shared annotated cell types found. Running DE on all cells combined.');
    shared_ct = "all_cells";
    ct1(:) = "all_cells";
    ct2(:) = "all_cells";
end

fprintf('Shared cell types (%d): %s\n', numel(shared_ct), strjoin(shared_ct, ', '));

% ---- Prepare output directory ---------------------------------------
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

% ---- DE per cell type -----------------------------------------------
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

    fprintf('DE for "%s": %d vs %d cells ... ', ct, n1, n2);

    T = [];
    try
        X = full(X1_all(:, mask1));
        Y = full(X2_all(:, mask2));
        T = sc_deg(X, Y, common_genes, 1, false);
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        continue;
    end

    % Post-process: label columns, handle Inf, build note table
    [T, Tnt] = pkg.in_DETableProcess(T, sample_id1, sample_id2, n1, n2);

    % Split into up/down regulated tables using standard thresholds
    [Tup, Tdn] = pkg.e_processdetable(T, paramset);

    fprintf('done. Up: %d  Down: %d  (FDR<0.01, |log2FC|≥1)\n', ...
        height(Tup), height(Tdn));

    ri = ri + 1;
    results(ri).cell_type = ct;
    results(ri).n1        = n1;
    results(ri).n2        = n2;
    results(ri).T         = T;
    results(ri).Tup       = Tup;
    results(ri).Tdn       = Tdn;

    % Save Excel file if out_dir provided
    if ~isempty(out_dir)
        outfile = sprintf('DE_%s_vs_%s_%s.xlsx', ...
            matlab.lang.makeValidName(sample_id1), ...
            matlab.lang.makeValidName(sample_id2), ...
            matlab.lang.makeValidName(string(ct)));
        filesaved = fullfile(out_dir, outfile);
        try
            writetable(T,   filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_genes');
            writetable(Tup, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Up-regulated');
            writetable(Tdn, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Down-regulated');
            writetable(Tnt, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Note');
            fprintf('  Saved: %s\n', filesaved);
        catch ME
            warning('Could not save %s: %s', filesaved, ME.message);
        end
    end
end

fprintf('\nDE analysis complete: %d cell type(s) analysed.\n', numel(results));

% ---- Print JSON summary for agent consumption -----------------------
% Everything printed here is captured by the MATLAB MCP tool and returned
% to the Python agent as text. The JSON block lets the Result Interpreting
% Agent reason about top DEGs without needing to re-enter MATLAB.
i_print_json_summary(results, sample_id1, sample_id2);
end


% ---- Helper: print JSON summary of top DEGs -------------------------
function i_print_json_summary(results, sample_id1, sample_id2, top_n)
if nargin < 4, top_n = 20; end

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
    'sample1',     char(sample_id1), ...
    'sample2',     char(sample_id2), ...
    'cell_types',  {cell_types});

fprintf('\n%%JSON_SUMMARY_BEGIN%%\n%s\n%%JSON_SUMMARY_END%%\n', ...
    jsonencode(summary, 'PrettyPrint', true));
end

function rows = i_table_to_structs(T, n)
rows = {};
if isempty(T), return; end
n = min(n, height(T));
for i = 1:n
    rows{end+1} = struct( ...
        'gene',      char(T.gene(i)), ...
        'log2FC',    round(T.avg_log2FC(i), 3), ...
        'p_val_adj', T.p_val_adj(i)); %#ok<AGROW>
end
end

% ---- Helper: locate and load cleandata.mat --------------------------
function sce = i_load_sce(sample_id, data_dir)
% Search data_dir/<any_study_id>/<sample_id>/cleandata.mat
hits = dir(fullfile(data_dir, '*', sample_id, 'cleandata.mat'));
if isempty(hits)
    % Flat fallback: data_dir/<sample_id>/cleandata.mat
    flat = fullfile(data_dir, sample_id, 'cleandata.mat');
    if isfile(flat)
        mat_path = flat;
    else
        error('llm:run_de_analysis:fileNotFound', ...
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
