function sce2 = toSCE2(obj)
%TOSCE2 Convert SingleCellExperiment to SingleCellExperiment2.
%
%   sce2 = toSCE2(sce)

    % Core matrix and gene names; cell IDs seeded from c_cell_id
    sce2 = SingleCellExperiment2(obj.X, obj.g, obj.c_cell_id);

    % --- Cell annotations -------------------------------------------------
    sce2.cellAnn.batch_id   = obj.c_batch_id(:);
    sce2.cellAnn.cluster_id = obj.c_cluster_id(:);
    sce2.cellAnn.cell_type  = obj.c_cell_type_tx(:);
    sce2.cellAnn.cell_cycle = obj.c_cell_cycle_tx(:);
    sce2.cellAnn.c          = obj.c(:);

    % Extra cell attributes from list_cell_attributes
    cnames = obj.list_cell_attributes(1:2:end);
    for k = 1:numel(cnames)
        fname = matlab.lang.makeValidName(cnames{k});
        sce2.cellAnn.(fname) = obj.list_cell_attributes{2*k}(:);
    end

    % --- Gene annotations -------------------------------------------------
    % geneAnn.gene_name already set by constructor
    gnames = obj.list_gene_attributes(1:2:end);
    for k = 1:numel(gnames)
        fname = matlab.lang.makeValidName(gnames{k});
        sce2.geneAnn.(fname) = obj.list_gene_attributes{2*k}(:);
    end

    % --- Embeddings -------------------------------------------------------
    sce2.embeddings = obj.struct_cell_embeddings;
    if ~isempty(obj.s)
        sce2.embeddings.s = obj.s;   % active / current embedding
    end

    % --- Clusterings ------------------------------------------------------
    sce2.clusterings = obj.struct_cell_clusterings;

    % --- Metadata ---------------------------------------------------------
    sce2.metadata.notes = obj.metadata;
end
