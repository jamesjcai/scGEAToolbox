function run_enrichr(results, out_dir, sample_id1, sample_id2, top_n, prefix)
% LLM.RUN_ENRICHR  Pathway enrichment (Enrichr) on DE results per cell type.
%
%   llm.run_enrichr(results, out_dir, sample_id1, sample_id2)
%   llm.run_enrichr(results, out_dir, sample_id1, sample_id2, top_n)
%
%   Runs Enrichr on the top up- and down-regulated genes for each cell type
%   in the results struct produced by llm.run_de_analysis. Uses the full
%   common gene list (results.T.gene) as the background.
%
%   Gene set libraries queried:
%     GO_Biological_Process_2025, GO_Molecular_Function_2025,
%     KEGG_2021_Human, Reactome_Pathways_2024
%
%   If out_dir and sample IDs are provided, enrichment tables are appended
%   as new sheets to the existing DE Excel files
%   (DE_<id1>_vs_<id2>_<celltype>.xlsx). If the file does not exist, a new
%   one is created containing only the Enrichr sheets.
%
%   A JSON summary of the top enriched terms per library is printed to
%   stdout between %ENRICHR_SUMMARY_BEGIN% and %ENRICHR_SUMMARY_END%
%   markers so the MATLAB MCP tool returns it to the agent.
%
%   Inputs:
%     results    - struct array from llm.run_de_analysis
%     out_dir    - folder for Excel output (pass [] to skip file writing)
%     sample_id1 - GSM accession of sample 1 (needed for Excel filename)
%     sample_id2 - GSM accession of sample 2 (needed for Excel filename)
%     top_n      - max genes submitted to Enrichr per direction (default 100)
%
%   Example (agent calls after run_de_analysis):
%     results = llm.run_de_analysis('GSM001', 'GSM002', data_dir, out_dir);
%     llm.run_enrichr(results, out_dir, 'GSM001', 'GSM002');

if nargin < 6 || isempty(prefix),     prefix     = "DE"; end
if nargin < 5 || isempty(top_n),      top_n      = 100; end
if nargin < 4 || isempty(sample_id2), sample_id2 = '';  end
if nargin < 3 || isempty(sample_id1), sample_id1 = '';  end
if nargin < 2,                         out_dir    = [];  end

genesets = ["GO_Biological_Process_2025", ...
            "GO_Molecular_Function_2025", ...
            "KEGG_2021_Human", ...
            "Reactome_Pathways_2024"];

lib_labels = ["GO_BP", "GO_MF", "KEGG", "Reactome"];

if ~isempty(out_dir) && ~isfolder(out_dir)
    mkdir(out_dir);
end

summary_ct = {};
errs = strings(0, 1);

for k = 1:numel(results)
    r = results(k);
    ct = r.cell_type;

    fprintf('\nEnrichr for "%s" ...', ct);

    if isempty(r.Tup) && isempty(r.Tdn)
        fprintf(' no DEGs, skipped.\n');
        continue;
    end

    background = r.T.gene;
    genes_up   = r.Tup.gene(1:min(top_n, height(r.Tup)));
    genes_dn   = r.Tdn.gene(1:min(top_n, height(r.Tdn)));

    % ---- Run Enrichr ------------------------------------------------
    [Tlist_up, errUp] = i_enrichrrun(genes_up, background, genesets, 'llm:run_enrichr:apiFailed');
    [Tlist_dn, errDn] = i_enrichrrun(genes_dn, background, genesets, 'llm:run_enrichr:apiFailed');
    errs = [errs; string(errUp); string(errDn)]; %#ok<AGROW>
    fprintf(' done.\n');

    % ---- Save to Excel ----------------------------------------------
    if ~isempty(out_dir) && ~isempty(sample_id1) && ~isempty(sample_id2)
        xlsfile = fullfile(out_dir, sprintf('%s_%s_vs_%s_%s.xlsx', ...
            prefix, ...
            matlab.lang.makeValidName(string(sample_id1)), ...
            matlab.lang.makeValidName(string(sample_id2)), ...
            matlab.lang.makeValidName(string(ct))));
        i_write_sheets(Tlist_up, xlsfile, lib_labels, 'Up');
        i_write_sheets(Tlist_dn, xlsfile, lib_labels, 'Dn');
        fprintf('  Appended Enrichr sheets to: %s\n', xlsfile);
    end

    % ---- Build summary struct for JSON output -----------------------
    up_struct = i_enrichrsummary(Tlist_up, lib_labels);
    dn_struct = i_enrichrsummary(Tlist_dn, lib_labels);

    summary_ct{end+1} = struct( ...
        'cell_type',  char(ct), ...
        'n_genes_up', numel(genes_up), ...
        'n_genes_dn', numel(genes_dn), ...
        'up',         up_struct, ...
        'dn',         dn_struct, ...
        'enrichr_error', i_joinErrors({errUp, errDn})); %#ok<AGROW>
end

% ---- Print JSON summary for agent -----------------------------------
% STATUS: 'completed', 'partial' (some Enrichr calls failed) or 'failed'
% (all did), each cell type carrying ENRICHR_ERROR, so an API outage is not
% read as "nothing enriched".
summary = struct('status', i_enrichrstatus(errs), 'cell_types', {summary_ct});
summaryJson = jsonencode(summary, 'PrettyPrint', true);
fprintf('\n%%ENRICHR_SUMMARY_BEGIN%%\n%s\n%%ENRICHR_SUMMARY_END%%\n', summaryJson);

% ---- Save JSON file so data_analysis_agent can read it back ---------
if ~isempty(out_dir) && ~isempty(sample_id1) && ~isempty(sample_id2)
    jsonFile = fullfile(out_dir, sprintf('%s_%s_vs_%s_enrichr.json', ...
        prefix, ...
        matlab.lang.makeValidName(string(sample_id1)), ...
        matlab.lang.makeValidName(string(sample_id2))));
    fid = fopen(jsonFile, 'w');
    if fid ~= -1
        fprintf(fid, '%s\n', summaryJson);
        fclose(fid);
    end
end

end


function msg = i_joinErrors(errs)
% The distinct non-empty messages, '' when there are none.
e = string(errs);
e = unique(e(e ~= ""));
msg = char(strjoin(e, "; "));
end


% ---- Append Enrichr tables as new Excel sheets ----------------------
function i_write_sheets(Tlist, xlsfile, lib_labels, direction)
for i = 1:numel(Tlist)
    T = Tlist{i};
    if isempty(T) || ~istable(T) || height(T) == 0, continue; end
    sheet = sprintf('%s_250_%s', direction, lib_labels(i));
    try
        writetable(T, xlsfile, 'FileType', 'spreadsheet', 'Sheet', sheet);
    catch ME
        warning('llm:run_enrichr:writeFailed', ...
            'Could not write sheet %s: %s', sheet, ME.message);
    end
end
end
