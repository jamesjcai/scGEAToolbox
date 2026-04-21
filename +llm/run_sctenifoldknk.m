function results = run_sctenifoldknk(sample_id, kogene, data_dir, out_dir)
% LLM.RUN_SCTENIFOLDKNK  Cell type-specific virtual gene knockout analysis.
%
%   results = llm.run_sctenifoldknk(sample_id, kogene)
%   results = llm.run_sctenifoldknk(sample_id, kogene, data_dir)
%   results = llm.run_sctenifoldknk(sample_id, kogene, data_dir, out_dir)
%
%   Builds a PCR-based Gene Regulatory Network (GRN) from one WT GEO
%   sample and performs a virtual knockout of the specified gene per cell
%   type. Uses manifold alignment to compare the WT GRN against the
%   pseudo-KO GRN to identify differentially regulated (DR) genes.
%
%   Method (lite — no tensor decomposition or subsampling):
%     1. Normalize and log-transform expression -> net.pcrnet -> GRN (A0)
%     2. Virtual KO: zero the row of kogene in A0 -> pseudo-KO GRN
%     3. Manifold alignment of A0 vs pseudo-KO GRN -> DR genes by
%        chi-squared test, FDR correction (BH). Genes with pAdjusted <
%        0.05 are considered differentially regulated by the KO.
%
%   Only one sample is required (WT). No experimental KO data is needed.
%
%   Inputs:
%     sample_id  - GSM accession of the WT sample (e.g. 'GSM2333580')
%     kogene     - gene to knock out virtually (string, e.g. 'Foxp3')
%     data_dir   - root folder containing downloaded .mat files, organised
%                  as data_dir/<study_id>/<sample_id>/cleandata.mat
%                  Default: 'data' (relative to MATLAB working directory).
%                  Pass an absolute path when calling from the agent.
%     out_dir    - folder to write Excel result files (optional).
%                  If omitted or empty, no files are written.
%
%   Output:
%     results - struct array with one element per analysed cell type:
%       .cell_type  - cell type label (string)
%       .n          - number of cells in this cell type
%       .kogene     - the knocked-out gene name
%       .T          - full DR table (all genes, sorted by drdist descending)
%                     Columns: sortid, genelist, drdist, FC, pValues, pAdjusted
%       .Tdr        - significant DR genes (pAdjusted < 0.05)
%
%   Reference:
%     Osorio D et al. (2022). scTenifoldKnk: An efficient virtual knockout
%     tool for gene function predictions via single-cell gene regulatory
%     network perturbation. Patterns 3(5):100434.
%     https://doi.org/10.1016/j.patter.2022.100434
%
%   Example (from agent via evaluate_matlab_code):
%     results = llm.run_sctenifoldknk('GSM2333580', 'Foxp3', ...
%                   'C:/abs/path/to/data', 'C:/abs/path/to/output');

if nargin < 3 || isempty(data_dir), data_dir = 'data'; end
if nargin < 4, out_dir = []; end

kogene = string(kogene);

results = struct('cell_type', {}, 'n', {}, 'kogene', {}, 'T', {}, 'Tdr', {});

% ---- Load SCE object ------------------------------------------------
sce = i_load_sce(sample_id, data_dir);
fprintf('Sample (%s): %d genes x %d cells\n', sample_id, sce.NumGenes, sce.NumCells);

% ---- Verify kogene exists in gene list ------------------------------
if ~any(sce.g == kogene)
    error('llm:run_sctenifoldknk:geneNotFound', ...
        'Gene "%s" not found in sample %s.', kogene, sample_id);
end
fprintf('Knockout gene: %s\n', kogene);

% ---- Limit to top 3000 HVGs for speed (preserve kogene) -------------
max_genes = 3000;
if sce.NumGenes > max_genes
    fprintf('Limiting to top %d HVGs (from %d genes)...\n', max_genes, sce.NumGenes);
    sce = i_limit_genes(sce, max_genes, kogene);
    fprintf('Genes after HVG selection: %d\n', sce.NumGenes);
end

% ---- QC filter ------------------------------------------------------
fprintf('Applying QC filter...\n');
sce = sce.qcfilterwhitelist(1000, 0.15, 15, 500, kogene);
fprintf('After QC: %d genes x %d cells\n', sce.NumGenes, sce.NumCells);

genelist = sce.g;
idx = find(genelist == kogene, 1);

% ---- Identify cell types --------------------------------------------
ct = sce.c_cell_type_tx;
cell_types = unique(ct);
cell_types = cell_types(~strcmpi(cell_types, 'undetermined'));

if isempty(cell_types)
    warning('llm:run_sctenifoldknk:noCellTypeAnnotation', ...
        'No annotated cell types found. Running on all cells combined.');
    cell_types = "all_cells";
    ct(:) = "all_cells";
end

fprintf('Cell types (%d): %s\n', numel(cell_types), strjoin(cell_types, ', '));

% ---- Prepare output directory ---------------------------------------
if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

i_log(out_dir, sprintf('START scTenifoldKnk (lite): %s | kogene=%s | %d cell types', ...
    sample_id, kogene, numel(cell_types)));

% ---- scTenifoldKnk (lite) per cell type -----------------------------
ri = 0;
for k = 1:numel(cell_types)
    ct_k = cell_types(k);
    mask = ct == ct_k;
    n = sum(mask);

    if n < 100
        fprintf('Skipping "%s": fewer than 100 cells (%d).\n', ct_k, n);
        i_log(out_dir, sprintf('SKIP "%s": n=%d (<100)', ct_k, n));
        continue;
    end

    fprintf('\nscTenifoldKnk (lite) for "%s": %d cells, KO=%s ...\n', ct_k, n, kogene);
    i_log(out_dir, sprintf('BEGIN "%s": n=%d kogene=%s', ct_k, n, kogene));

    X = full(sce.X(:, mask));
    X = sc_norm(X);
    X = log1p(X);

    useGPU = pkg.i_usegpu(X);
    useparallel = ~useGPU;

    A0 = [];
    try
        disp('Constructing gene regulatory network...')
        A0 = net.pcrnet(X, 3, true, true, useparallel, ~useparallel, useGPU);
    catch ME
        fprintf('FAILED (network construction): %s\n', ME.message);
        i_log(out_dir, sprintf('FAILED "%s" (pcrnet): %s', ct_k, ME.message));
        continue;
    end

    nlinks = nnz(A0(idx, :) ~= 0);
    if nlinks == 0
        fprintf('Skipping "%s": KO gene (%s) has no links in network.\n', ct_k, kogene);
        i_log(out_dir, sprintf('SKIP "%s": %s has no links', ct_k, kogene));
        continue;
    end
    if nlinks < 50
        fprintf('Note: KO gene (%s) has few links (n=%d) in "%s".\n', kogene, nlinks, ct_k);
    end

    T = [];
    try
        disp('>> [T] = ten.i_knk(A0, idx, genelist, true);')
        T = ten.i_knk(A0, idx, genelist, true);
    catch ME
        fprintf('FAILED (knockout): %s\n', ME.message);
        i_log(out_dir, sprintf('FAILED "%s" (i_knk): %s', ct_k, ME.message));
        continue;
    end

    if isempty(T) || ~istable(T) || ~ismember('pAdjusted', T.Properties.VariableNames)
        fprintf('  No valid results returned.\n');
        i_log(out_dir, sprintf('FAILED "%s": no valid results returned', ct_k));
        continue;
    end

    T = sortrows(T, 'drdist', 'descend');
    Tdr = T(T.pAdjusted < 0.05, :);
    fprintf('  DR genes (pAdj<0.05): %d\n', height(Tdr));
    i_log(out_dir, sprintf('DONE "%s": %d DR genes (pAdj<0.05)', ct_k, height(Tdr)));

    ri = ri + 1;
    results(ri).cell_type = ct_k;
    results(ri).n         = n;
    results(ri).kogene    = kogene;
    results(ri).T         = T;
    results(ri).Tdr       = Tdr;

    if ~isempty(out_dir)
        outfile = sprintf('KNK_%s_%s_%s.xlsx', ...
            matlab.lang.makeValidName(sample_id), ...
            matlab.lang.makeValidName(kogene), ...
            matlab.lang.makeValidName(string(ct_k)));
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

fprintf('\nscTenifoldKnk complete: %d cell type(s) analysed.\n', numel(results));
i_log(out_dir, sprintf('COMPLETE scTenifoldKnk: %d cell type(s) analysed; emitting JSON summary', numel(results)));

% ---- Print JSON summary for agent consumption -----------------------
i_print_json_summary(results, sample_id, kogene);
end


% ---- Helper: print JSON summary of top DR genes ---------------------
function i_print_json_summary(results, sample_id, kogene, top_n)
if nargin < 4, top_n = 20; end

cell_types = {};
for k = 1:numel(results)
    r = results(k);
    top_dr = i_table_to_structs(r.Tdr, top_n);
    cell_types{end+1} = struct( ...
        'cell_type',  char(r.cell_type), ...
        'n',          r.n, ...
        'n_dr_genes', height(r.Tdr), ...
        'top_dr',     {top_dr}); %#ok<AGROW>
end

summary = struct( ...
    'sample',     char(sample_id), ...
    'kogene',     char(kogene), ...
    'cell_types', {cell_types});

fprintf('\n%%SCTENIFOLDKNK_JSON_SUMMARY_BEGIN%%\n%s\n%%SCTENIFOLDKNK_JSON_SUMMARY_END%%\n', ...
    jsonencode(summary, 'PrettyPrint', true));
end


% ---- Helper: extract top rows as struct list for JSON ---------------
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
hits = dir(fullfile(data_dir, '*', sample_id, 'cleandata.mat'));
if isempty(hits)
    flat = fullfile(data_dir, sample_id, 'cleandata.mat');
    if isfile(flat)
        mat_path = flat;
    else
        error('llm:run_sctenifoldknk:fileNotFound', ...
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
end
end


% ---- Helper: keep top n HVGs; always include mustkeep gene ----------
function sce = i_limit_genes(sce, n, mustkeep)
T_hvg = sc_splinefit(sce.X, sce.g);
glist = T_hvg.genes(1:min(n, sce.NumGenes));
if nargin >= 3 && ~isempty(mustkeep) && ~ismember(mustkeep, glist)
    glist(end) = mustkeep;
end
[~, hidx] = ismember(glist, sce.g);
hidx = hidx(hidx > 0);
sce.X = sce.X(hidx, :);
sce.g = sce.g(hidx);
end
