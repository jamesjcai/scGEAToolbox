function results = run_sctenifoldnet(sample_id1, sample_id2, data_dir, out_dir)
% LLM.RUN_SCTENIFOLDNET  Cell type-specific scTenifoldNet GRN comparison.
%
%   results = llm.run_sctenifoldnet(sample_id1, sample_id2)
%   results = llm.run_sctenifoldnet(sample_id1, sample_id2, data_dir)
%   results = llm.run_sctenifoldnet(sample_id1, sample_id2, data_dir, out_dir)
%
%   Builds PCR-based Gene Regulatory Networks (GRNs) for two GEO samples
%   and identifies differentially regulated (DR) genes per cell type.
%
%   Method: ten.sctenifoldnet with its defaults, once per shared cell type, on
%   raw counts. That is, per sample:
%     1. Drop ribosomal genes; libsize-normalize and log1p
%     2. Build 10 PCNets, each on 500 bootstrap-resampled cells (ten.i_nc)
%     3. Denoise them into one network by CP tensor decomposition
%   then symmetrize both networks, align them (ten.i_ma), and call DR genes
%   with a chi-squared test on squared alignment distances, BH-corrected
%   (ten.i_dr). Genes with pAdjusted < 0.05 are DR genes.
%
%   Requires the Tensor Toolbox (Setup > Install Tensor Toolbox); without it
%   every cell type fails with the toolbox's install message.
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
%       .T          - full DR table (all genes, sorted by drdist descending)
%                     Columns: sortid, genelist, drdist, FC, pValues, pAdjusted
%       .Tdr        - significant DR genes (pAdjusted < 0.05)
%
%   Reference:
%     Osorio D et al. (2022). scTenifoldNet: A Machine Learning Workflow for
%     Constructing and Comparing Transcriptome-wide Gene Regulatory Networks
%     from Single-Cell Data. Patterns 3(3):100434.
%     https://doi.org/10.1016/j.patter.2020.100434
%
%   Example (from agent via evaluate_matlab_code):
%     results = llm.run_sctenifoldnet('GSM2333580', 'GSM2333581', ...
%                   'C:/abs/path/to/data', 'C:/abs/path/to/output');

if nargin < 3 || isempty(data_dir), data_dir = 'data'; end
if nargin < 4, out_dir = []; end

results = struct('cell_type', {}, 'n1', {}, 'n2', {}, 'T', {}, 'Tdr', {});

% ---- Load SCE objects -----------------------------------------------
sce1 = i_load_sce(sample_id1, data_dir);
sce2 = i_load_sce(sample_id2, data_dir);

fprintf('Sample 1 (%s): %d genes x %d cells\n', sample_id1, sce1.NumGenes, sce1.NumCells);
fprintf('Sample 2 (%s): %d genes x %d cells\n', sample_id2, sce2.NumGenes, sce2.NumCells);

% ---- Align to common gene set ---------------------------------------
[common_genes, idx1, idx2] = intersect(sce1.g, sce2.g, 'stable');
if isempty(common_genes)
    error('llm:run_sctenifoldnet:noCommonGenes', ...
        'No common genes between %s and %s.', sample_id1, sample_id2);
end
fprintf('Common genes: %d\n', numel(common_genes));

X1_all = sce1.X(idx1, :);
X2_all = sce2.X(idx2, :);

% ---- Limit to top 3000 HVGs for speed -------------------------------
max_genes = 3000;
if numel(common_genes) > max_genes
    fprintf('Limiting to top %d HVGs (from %d common genes)...\n', max_genes, numel(common_genes));
    sce_tmp = SingleCellExperiment(X1_all, common_genes);
    sce_tmp = llm.i_limit_genes(sce_tmp, max_genes);
    [~, hidx] = ismember(sce_tmp.g, common_genes);
    common_genes = sce_tmp.g;
    X1_all = X1_all(hidx, :);
    X2_all = X2_all(hidx, :);
    fprintf('Genes after HVG selection: %d\n', numel(common_genes));
end

% ---- QC filter on the gene-limited data -----------------------------
sce1.X = X1_all; sce1.g = common_genes;
sce2.X = X2_all; sce2.g = common_genes;
fprintf('Applying QC filter (sample 1)...\n');
sce1 = sce1.qcfilterwhitelist(1000, 0.15, 15, 500, []);
fprintf('Applying QC filter (sample 2)...\n');
sce2 = sce2.qcfilterwhitelist(1000, 0.15, 15, 500, []);
% QCFILTERWHITELIST drops genes detected in fewer than 15 cells, sample by
% sample, so the two gene lists can differ afterwards. They were used as if
% they still matched COMMON_GENES: unequal lengths made every cell type fail
% in i_ma, and equal lengths with different genes misaligned the rows and
% named the wrong genes as differentially regulated.
[common_genes, i1, i2] = intersect(sce1.g, sce2.g, 'stable');
X1_all = sce1.X(i1, :);
X2_all = sce2.X(i2, :);
fprintf('Genes kept by QC in both samples: %d\n', numel(common_genes));
fprintf('After QC: %s (sample 1), %s (sample 2)\n', ...
    pkg.i_plural(sce1.NumCells, 'cell'), ...
    pkg.i_plural(sce2.NumCells, 'cell'));

% ---- Identify shared cell types -------------------------------------
ct1 = sce1.c_cell_type_tx;
ct2 = sce2.c_cell_type_tx;

shared_ct = intersect(unique(ct1), unique(ct2));
shared_ct = shared_ct(~strcmpi(shared_ct, 'undetermined'));

if isempty(shared_ct)
    warning('llm:run_sctenifoldnet:noCellTypeAnnotation', ...
        'No shared annotated cell types found. Running on all cells combined.');
    shared_ct = "all_cells";
    ct1(:) = "all_cells";
    ct2(:) = "all_cells";
end

fprintf('Shared cell types (%d): %s\n', numel(shared_ct), strjoin(shared_ct, ', '));

% ---- Prepare output directory ---------------------------------------
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

i_log(out_dir, sprintf('START scTenifoldNet: %s vs %s | %d cell types', ...
    sample_id1, sample_id2, numel(shared_ct)));

% ---- scTenifoldNet per cell type ------------------------------------
ri = 0;
nTooSmall = 0;
nFailed = 0;
for k = 1:numel(shared_ct)
    ct = shared_ct(k);
    mask1 = ct1 == ct;
    mask2 = ct2 == ct;
    n1 = sum(mask1);
    n2 = sum(mask2);

    if n1 < 100 || n2 < 100
        fprintf('Skipping "%s": fewer than 100 cells (%d in sample1, %d in sample2).\n', ...
            ct, n1, n2);
        i_log(out_dir, sprintf('SKIP "%s": n1=%d n2=%d (both need >=100)', ct, n1, n2));
        nTooSmall = nTooSmall + 1;
        continue;
    end

    fprintf('\nscTenifoldNet for "%s": %d vs %d cells ...\n', ct, n1, n2);
    i_log(out_dir, sprintf('BEGIN "%s": n1=%d n2=%d', ct, n1, n2));

    % All cells go in: ten.sctenifoldnet draws its own 500-cell bootstrap
    % subsamples per network, so a cap here would only discard cells. It
    % normalizes the raw counts itself, and savegrn=false keeps it from
    % writing A0_*.mat/A1_*.mat into the working folder.
    try
        T = ten.sctenifoldnet(full(X1_all(:, mask1)), full(X2_all(:, mask2)), ...
            common_genes, 'savegrn', false);
    catch ME
        fprintf('FAILED: %s\n', ME.message);
        i_log(out_dir, sprintf('FAILED "%s": %s', ct, ME.message));
        nFailed = nFailed + 1;
        continue;
    end

    if isempty(T) || ~istable(T) || ~ismember('pAdjusted', T.Properties.VariableNames)
        fprintf('  No valid results returned.\n');
        i_log(out_dir, sprintf('FAILED "%s": no valid results returned', ct));
        nFailed = nFailed + 1;
        continue;
    end

    T = sortrows(T, 'drdist', 'descend');
    Tdr = T(T.pAdjusted < 0.05, :);
    fprintf('  DR genes (pAdj<0.05): %d\n', height(Tdr));
    i_log(out_dir, sprintf('DONE "%s": %d DR genes (pAdj<0.05)', ct, height(Tdr)));

    ri = ri + 1;
    results(ri).cell_type = ct;
    results(ri).n1        = n1;
    results(ri).n2        = n2;
    results(ri).T         = T;
    results(ri).Tdr       = Tdr;

    if ~isempty(out_dir)
        outfile = sprintf('scTenifoldNet_%s_vs_%s_%s.xlsx', ...
            matlab.lang.makeValidName(sample_id1), ...
            matlab.lang.makeValidName(sample_id2), ...
            matlab.lang.makeValidName(string(ct)));
        filesaved = fullfile(out_dir, outfile);
        try
            writetable(Tdr, filesaved, 'FileType', 'spreadsheet', 'Sheet', 'DR_genes');
            writetable(T,   filesaved, 'FileType', 'spreadsheet', 'Sheet', 'All_genes');
            fprintf('  Saved: %s\n', filesaved);
        catch ME
            warning('Could not save %s: %s', filesaved, ME.message);
        end
    end
end

fprintf('\nscTenifoldNet complete: %d cell type(s) analysed.\n', numel(results));
i_log(out_dir, sprintf('COMPLETE scTenifoldNet: %d cell type(s) analysed; emitting JSON summary', numel(results)));

% ---- Print JSON summary for agent consumption -----------------------
if isempty(results)
    % Say which: this always blamed the 100-cell minimum, including when
    % every cell type had failed in the analysis itself.
    no_results_reason = sprintf(['No cell type produced results: %d had ' ...
        'fewer than 100 cells in a sample, %d failed in the analysis (see ' ...
        'the run log).'], nTooSmall, nFailed);
else
    no_results_reason = '';
end
i_print_json_summary(results, sample_id1, sample_id2, no_results_reason);
end


% ---- Helper: print JSON summary of top DR genes ---------------------
function i_print_json_summary(results, sample_id1, sample_id2, no_results_reason, top_n)
if nargin < 5, top_n = 20; end
if nargin < 4, no_results_reason = ''; end

% top_dr is taken from the FULL ranked T table (sorted by drdist desc),
% NOT from Tdr (pAdj<0.05 subset), so it matches the gene list used for
% Enrichr pathway enrichment. n_dr_significant reports the significance count.
cell_types = {};
for k = 1:numel(results)
    r = results(k);
    top_dr = i_table_to_structs(r.T, top_n);
    cell_types{end+1} = struct( ...
        'cell_type',        char(r.cell_type), ...
        'n1',               r.n1, ...
        'n2',               r.n2, ...
        'n_dr_significant', height(r.Tdr), ...
        'n_dr_total',       height(r.T), ...
        'top_dr',           {top_dr}); %#ok<AGROW>
end

if isempty(cell_types)
    reason = no_results_reason;
    if isempty(reason)
        reason = 'scTenifoldNet produced no results.';
    end
    summary = struct( ...
        'sample1',    char(sample_id1), ...
        'sample2',    char(sample_id2), ...
        'status',     'no_results', ...
        'reason',     reason, ...
        'cell_types', {cell_types});
else
    summary = struct( ...
        'sample1',    char(sample_id1), ...
        'sample2',    char(sample_id2), ...
        'status',     'completed', ...
        'method',     'scTenifoldNet (ten.sctenifoldnet defaults: 10 bootstrap PCNets per sample, CP tensor denoising, manifold alignment)', ...
        'cell_types', {cell_types});
end

fprintf('\n%%SCTENIFOLDNET_JSON_SUMMARY_BEGIN%%\n%s\n%%SCTENIFOLDNET_JSON_SUMMARY_END%%\n', ...
    jsonencode(summary, 'PrettyPrint', true));
end


% ---- Helper: extract top rows as struct list for JSON ---------------
% NOTE: FC here is network fold change (ratio of regulatory weights in the
% aligned manifold), NOT expression fold change. drdist is the Euclidean
% distance between a gene's positions in the two aligned manifolds.
function rows = i_table_to_structs(T, n)
rows = {};
if isempty(T), return; end
n = min(n, height(T));
for i = 1:n
    rows{end+1} = struct( ...
        'gene',      char(T.genelist(i)), ...
        'drdist',    round(T.drdist(i), 6), ...
        'FC',        round(T.FC(i), 4), ...
        'pAdjusted', T.pAdjusted(i)); %#ok<AGROW>
end
end


% ---- Helper: locate and load cleandata.mat --------------------------
function sce = i_load_sce(sample_id, data_dir)
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


