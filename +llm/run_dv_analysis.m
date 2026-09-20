function results = run_dv_analysis(sample_id1, sample_id2, data_dir, out_dir, method, direction)
% LLM.RUN_DV_ANALYSIS  Cell type-specific DV analysis between two GEO samples.
%
%   results = llm.run_dv_analysis(sample_id1, sample_id2)
%   results = llm.run_dv_analysis(sample_id1, sample_id2, data_dir)
%   results = llm.run_dv_analysis(sample_id1, sample_id2, data_dir, out_dir)
%   results = llm.run_dv_analysis(..., out_dir, method)
%   results = llm.run_dv_analysis(..., out_dir, method, direction)
%
%   Loads cleandata.mat for each sample (each file contains a
%   SingleCellExperiment variable named 'sce') and performs differential
%   variability (DV) analysis with SC_DVG for every cell type shared
%   between the two samples.
%
%   DiffSign interpretation (set by DIRECTION; the JSON summary names it):
%     'mean' (default)
%       > 0  — up-regulated: mean expression higher in sample 1 (tested)
%       < 0  — down-regulated: mean expression lower in sample 1 than in
%              sample 2 (baseline)
%     'deviation'
%       > 0  — higher transcriptional variability in sample 1
%       < 0  — lower transcriptional variability in sample 1
%   DiffDist, not DiffSign, carries the size of the variability difference.
%
%   If out_dir is provided, results are saved as Excel files
%   (DV_<id1>_vs_<id2>_<celltype>.xlsx; a non-default method adds
%   _<method> after DV, and direction 'deviation' adds _devsign) with
%   sheets: All genes, Up-regulated, Down-regulated, Note ('Variability
%   increasing'/'decreasing' in place of up/down for 'deviation').
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
%     method     - reference curve each sample is scored against:
%                  'splinefit' (default) fits a smoothing spline;
%                  'analytic' uses the closed-form gamma-Poisson curve,
%                  which is defined at every mean and so discards no genes
%                  at the ends of the fitted range. Both return the same
%                  table, so the rest of this function is unaffected.
%     direction  - how DV genes are split into up and down (see SC_DVG):
%                  'mean' (default) by mean expression, or 'deviation'
%                  by deviation from each sample's curve.
%
%   Output:
%     results - struct array with one element per shared cell type:
%       .cell_type  - cell type label (string)
%       .n1         - number of cells from sample 1
%       .n2         - number of cells from sample 2
%       .T          - full DV table (all genes, sorted by DiffDist descending)
%       .Tup        - significant DV genes up in sample 1 (DiffSign > 0)
%       .Tdn        - significant DV genes down in sample 1 (DiffSign < 0)
%
%   Example (from agent via evaluate_matlab_code):
%     results = llm.run_dv_analysis('GSM2333580', 'GSM2333581', ...
%                   'C:/abs/path/to/data', 'C:/abs/path/to/output');

if nargin < 3 || isempty(data_dir), data_dir = 'data'; end
if nargin < 4, out_dir = []; end
if nargin < 5 || isempty(method), method = 'splinefit'; end
method = validatestring(method, {'splinefit', 'analytic'}, ...
    mfilename, 'method', 5);
if nargin < 6 || isempty(direction), direction = 'mean'; end
direction = validatestring(direction, {'mean', 'deviation'}, ...
    mfilename, 'direction', 6);

max_cells = 2000;   % subsample per cell type for speed

results = struct('cell_type', {}, 'n1', {}, 'n2', {}, ...
                 'T', {}, 'Tup', {}, 'Tdn', {});
skipped = {};   % records of skipped cell types with reasons

% ---- Load SCE objects -----------------------------------------------
sce1 = i_load_sce(sample_id1, data_dir);
sce2 = i_load_sce(sample_id2, data_dir);

fprintf('Sample 1 (%s): %d genes x %d cells\n', sample_id1, sce1.NumGenes, sce1.NumCells);
fprintf('Sample 2 (%s): %d genes x %d cells\n', sample_id2, sce2.NumGenes, sce2.NumCells);

% ---- Identify shared cell types -------------------------------------
ct1 = sce1.c_cell_type_tx;
ct2 = sce2.c_cell_type_tx;

shared_ct = intersect(unique(ct1), unique(ct2));
shared_ct = shared_ct(~strcmpi(shared_ct, 'undetermined'));

if isempty(shared_ct)
    warning('llm:run_dv_analysis:noCellTypeAnnotation', ...
        'No shared annotated cell types found. Running DV on all cells combined.');
    shared_ct = "all_cells";
    ct1(:) = "all_cells";
    ct2(:) = "all_cells";
end

fprintf('Shared cell types (%d): %s\n', numel(shared_ct), strjoin(shared_ct, ', '));

% ---- Prepare output directory ---------------------------------------
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

% ---- DV per cell type -----------------------------------------------
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
        skipped{end+1} = struct('cell_type', char(ct), 'n1', n1, 'n2', n2, ...
            'reason', sprintf('insufficient_cells (need ≥500 each; have %d vs %d)', n1, n2));
        continue;
    end

    fprintf('DV for "%s": %d vs %d cells ... ', ct, n1, n2);

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

    % Subset SCE objects by cell type
    sce1_ct = sce1.selectcells(idx1);
    sce2_ct = sce2.selectcells(idx2);

    % QC filter (remove low-quality cells/genes after subsetting)
    sce1_ct = sce1_ct.qcfilter;
    sce2_ct = sce2_ct.qcfilter;

    if sce1_ct.NumCells < 10 || sce2_ct.NumCells < 10 || ...
       sce1_ct.NumGenes < 10 || sce2_ct.NumGenes < 10
        fprintf('SKIPPED (too few cells/genes after QC filter).\n');
        skipped{end+1} = struct('cell_type', char(ct), 'n1', n1, 'n2', n2, ...
            'reason', 'insufficient_cells_after_qc_filter');
        continue;
    end

    T = [];
    try
        T = sc_dvg(sce1_ct, sce2_ct, ...
            {char(sample_id1)}, {char(sample_id2)}, method, direction);
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        skipped{end+1} = struct('cell_type', char(ct), 'n1', n1, 'n2', n2, ...
            'reason', sprintf('analysis_error: %s', ME.message));
        continue;
    end

    % Label columns (adds sample-specific headers and note fields)
    [T, Tnt] = pkg.in_DVTableProcess(T, {char(sample_id1)}, {char(sample_id2)}, direction);

    % Split significant genes (pval < 0.05) into up (DiffSign > 0, sample 1) and down
    Tup = T(T.DiffSign > 0 & T.pval < 0.05, :);
    Tdn = T(T.DiffSign < 0 & T.pval < 0.05, :);

    fprintf('done. Up: %d  Down: %d  (pval<0.05)\n', height(Tup), height(Tdn));

    ri = ri + 1;
    results(ri).cell_type = ct;
    results(ri).n1        = n1;
    results(ri).n2        = n2;
    results(ri).T         = T;
    results(ri).Tup       = Tup;
    results(ri).Tdn       = Tdn;

    % Save Excel file if out_dir provided
    if ~isempty(out_dir)
        % A non-default curve is named in the file, so that a second
        % run does not overwrite the first. The default keeps the name
        % it has always written, which callers look for.
        methodinfix = '';
        if ~strcmp(method, 'splinefit')
            methodinfix = ['_', method];
        end
        uplabel = 'Up-regulated (p<0.05)';
        dnlabel = 'Down-regulated (p<0.05)';
        if strcmp(direction, 'deviation')
            methodinfix = [methodinfix, '_devsign'];
            % 31 characters, Excel's sheet-name limit.
            uplabel = 'Variability increasing (p<0.05)';
            dnlabel = 'Variability decreasing (p<0.05)';
        end
        outfile = sprintf('DV%s_%s_vs_%s_%s.xlsx', methodinfix, ...
            matlab.lang.makeValidName(sample_id1), ...
            matlab.lang.makeValidName(sample_id2), ...
            matlab.lang.makeValidName(string(ct)));
        filesaved = fullfile(out_dir, outfile);
        try
            writetable(T,   filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All genes');
            writetable(Tup, filesaved, 'FileType', 'spreadsheet', 'Sheet', uplabel);
            writetable(Tdn, filesaved, 'FileType', 'spreadsheet', 'Sheet', dnlabel);
            writetable(Tnt, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'Note');
            fprintf('  Saved: %s\n', filesaved);
        catch ME
            warning('Could not save %s: %s', filesaved, ME.message);
        end
    end
end

fprintf('\nDV analysis complete: %d cell type(s) analysed, %d skipped.\n', ...
    numel(results), numel(skipped));

% ---- Print JSON summary for agent consumption -----------------------
i_print_json_summary(results, skipped, sample_id1, sample_id2, method, direction);
end


% ---- Helper: print JSON summary of top DV genes ---------------------
function i_print_json_summary(results, skipped, sample_id1, sample_id2, ...
        method, direction, top_n)
if nargin < 7, top_n = 20; end

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

% Build no_results_reason when nothing passed
no_results_reason = '';
if isempty(cell_types) && ~isempty(skipped)
    no_results_reason = sprintf(['All %d shared cell type(s) were skipped. ' ...
        'DV requires ≥500 cells per cell type in each sample. See skipped[] ' ...
        'for per-cell-type counts and reasons.'], numel(skipped));
elseif isempty(cell_types)
    no_results_reason = 'No shared annotated cell types found between the two samples.';
end

summary = struct( ...
    'sample1',           char(sample_id1), ...
    'sample2',           char(sample_id2), ...
    'method',            char(method), ...
    'direction',         char(direction), ...
    'cell_types',        {cell_types}, ...
    'skipped',           {skipped}, ...
    'no_results_reason', no_results_reason);

fprintf('\n%%DV_JSON_SUMMARY_BEGIN%%\n%s\n%%DV_JSON_SUMMARY_END%%\n', ...
    jsonencode(summary, 'PrettyPrint', true));
end


% ---- Helper: extract top rows as struct list for JSON ---------------
function rows = i_table_to_structs(T, n)
rows = {};
if isempty(T), return; end
n = min(n, height(T));
for i = 1:n
    rows{end+1} = struct( ...
        'gene',      char(T.gene(i)), ...
        'DiffDist',  round(T.DiffDist(i), 4), ...
        'pval',      T.pval(i)); %#ok<AGROW>
end
end


% ---- Helper: locate and load cleandata.mat --------------------------
function sce = i_load_sce(sample_id, data_dir)
sce = llm.i_load_sce(sample_id, data_dir);
end
