function exportToJsonl(sce, outfile, dataset_id, titleText)
%EXPORTTOJSONL Export SingleCellExperiment to JSONL for vector DB
%
%  exportToJsonl(sce, 'out.jsonl', 'GSM12345', 'Breast cancer study');

if nargin < 4
    titleText = '';
end

% Basic stats
NumGenes = size(sce.X, 1);
NumCells = size(sce.X, 2);

% Gene list summary (first 20 genes or fewer)
gene_list_summary = sce.g(1:min(20, NumGenes));

% ---- Cell types ----
cell_types = struct();
if isprop(sce, 'c_cell_type_tx') && ~isempty(sce.c_cell_type_tx)
    [utypes, ~, idx] = unique(sce.c_cell_type_tx);
    counts = accumarray(idx, 1);
    cell_types.counts = containers.Map(utypes, counts);
    cell_types.source = 'c_cell_type_tx';
end

% ---- Batches ----
batches = struct();
if isprop(sce, 'c_batch_id') && ~isempty(sce.c_batch_id)
    [ubatch, ~, idx] = unique(sce.c_batch_id);
    counts = accumarray(idx, 1);
    batches.unique_ids = ubatch;
    batches.distribution = containers.Map(ubatch, counts);
    batches.source = 'c_batch_id';
end

% ---- Clusters ----
clusters = struct();
if isprop(sce, 'c_cluster_id') && ~isempty(sce.c_cluster_id)
    clusters.n_clusters = numel(unique(sce.c_cluster_id));
    clusters.source = 'c_cluster_id';
end

% ---- Cell cycle ----
attributes = struct();
if isprop(sce, 'c_cell_cycle_tx') && ~isempty(sce.c_cell_cycle_tx)
    [ucycle, ~, idx] = unique(sce.c_cell_cycle_tx);
    counts = accumarray(idx, 1);
    attributes.cell_cycle = containers.Map(ucycle, counts);
end

% ---- Metadata ----
metadata = struct();
if isprop(sce, 'metadata') && isstruct(sce.metadata)
    metadata = sce.metadata;
end
if isprop(sce, 'notes') && ~isempty(sce.notes)
    metadata.notes = sce.notes;
end
metadata.created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');

% ---- Embeddings (optional: UMAP/TSNE if available) ----
% embeddings = struct();
% if isprop(sce, 'struct_cell_embeddings') && ~isempty(sce.struct_cell_embeddings)
%     fns = fieldnames(sce.struct_cell_embeddings);
%     for i = 1:numel(fns)
%         mat = sce.struct_cell_embeddings.(fns{i});
%         if size(mat,2) <= 3   % only keep low-dim embeddings
%             embeddings.(fns{i}) = mat(1:min(50, size(mat,1)), :); % preview only
%         end
%     end
% end

% ---- Construct JSON struct ----
record = struct( ...
    'dataset_id', dataset_id, ...
    'title', titleText, ...
    'NumGenes', NumGenes, ...
    'NumCells', NumCells, ...
    'gene_list_summary', gene_list_summary, ...
    'cell_types', cell_types, ...
    'clusters', clusters, ...
    'batches', batches, ...
    'attributes', attributes, ...
    'metadata', metadata ...
    );

% ---- Encode and write ----
jsonstr = jsonencode(record, 'PrettyPrint', true);

fid = fopen(outfile, 'a'); % append mode
fprintf(fid, '%s\n', jsonstr);
fclose(fid);

fprintf('✅ Exported dataset %s to %s\n', dataset_id, outfile);

end
