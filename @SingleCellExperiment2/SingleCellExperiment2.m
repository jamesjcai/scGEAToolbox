classdef SingleCellExperiment2 < handle & matlab.mixin.Copyable
    % SingleCellExperiment2
    % A MATLAB class inspired by AnnData (.h5ad) structure
    
    properties
        % --- Core ---
        X                  % main data matrix (genes × cells)
        
        % --- Annotations ---
        cellAnn table      % like obs (per-cell annotations)
        geneAnn table      % like var (per-gene annotations)
        
        % --- Optional ---
        layers struct      % alternative matrices (raw, lognorm, scaled, ...)
        embeddings struct  % low-dim representations (PCA, UMAP, tSNE)
        clusterings struct % clustering results
        metadata struct    % unstructured metadata (colors, settings, etc.)
    end

    properties (Dependent, SetAccess=private)
        NumCells
        NumGenes
    end    
    
    methods
        %% Constructor
        function obj = SingleCellExperiment2(X, g, cellIDs)
            if nargin > 0
                % Core matrix
                obj.X = X;
                
                % Gene annotations (var)
                if nargin >= 2 && ~isempty(g)
                    obj.geneAnn = table(g(:), 'VariableNames', "gene_name");
                else
                    obj.geneAnn = table( ...
                        "Gene" + (1:size(X,1))', ...
                        'VariableNames', "gene_name");
                end
                
                % Cell annotations (obs)
                if nargin >= 3 && ~isempty(cellIDs)
                    obj.cellAnn = table(cellIDs(:), 'VariableNames', "cell_id");
                else
                    obj.cellAnn = table( ...
                        "Cell" + (1:size(X,2))', ...
                        'VariableNames', "cell_id");
                end
                
                % Empty slots
                obj.layers     = struct();
                obj.embeddings = struct();
                obj.clusterings= struct();
                obj.metadata   = struct();
            else
                % empty object
                obj.X = [];
                obj.geneAnn = table();
                obj.cellAnn = table();
                obj.layers = struct();
                obj.embeddings = struct();
                obj.clusterings = struct();
                obj.metadata = struct();
            end
        end

        %% ---- Dependent properties ----
        function n = get.NumCells(obj)
            n = size(obj.X,2);
        end
        function n = get.NumGenes(obj)
            n = size(obj.X,1);
        end

        %% --- Subset cells by indices or logical mask
        function obj = subsetCells(obj, idx)
            if islogical(idx)
                if numel(idx) ~= obj.NumCells
                    error("Logical mask length must equal number of cells");
                end
            elseif isnumeric(idx)
                if max(idx) > obj.NumCells
                    error("Cell index out of range");
                end
            else
                error("idx must be numeric or logical");
            end
    
            obj.X       = obj.X(:, idx);
            obj.cellAnn = obj.cellAnn(idx, :);
    
            % Subset embeddings
            fn = fieldnames(obj.embeddings);
            for k = 1:numel(fn)
                val = obj.embeddings.(fn{k});
                if size(val,1) == obj.NumCells
                    obj.embeddings.(fn{k}) = val(idx, :);
                end
            end
        end
    
        %% --- Subset genes by indices or logical mask
        function obj = subsetGenes(obj, idx)
            if islogical(idx)
                if numel(idx) ~= obj.NumGenes
                    error("Logical mask length must equal number of genes");
                end
            elseif isnumeric(idx)
                if max(idx) > obj.NumGenes
                    error("Gene index out of range");
                end
            else
                error("idx must be numeric or logical");
            end
    
            obj.X       = obj.X(idx, :);
            obj.geneAnn = structfun(@(x) x(idx,:), obj.geneAnn, 'UniformOutput', false);
        end
    
        %% --- Remove empty cells (all-zero expression)
        function obj = rmEmptyCells(obj)
            keep = sum(obj.X,1) > 0;
            obj  = obj.subsetCells(keep);
        end
    
        %% --- Remove empty genes (all-zero expression)
        function obj = rmEmptyGenes(obj)
            keep = sum(obj.X,2) > 0;
            obj  = obj.subsetGenes(keep);
        end
    
        %% --- Keep only highly expressed genes
        function obj = filterGenesByMinCells(obj, minCells)
            if nargin < 2, minCells = 10; end
            keep = sum(obj.X > 0, 2) >= minCells;
            obj  = obj.subsetGenes(keep);
        end
    
        %% --- Keep only cells with library size above cutoff
        function obj = filterCellsByLibSize(obj, minCounts)
            if nargin < 2, minCounts = 1000; end
            keep = sum(obj.X, 1) >= minCounts;
            obj  = obj.subsetCells(keep);
        end

        %% Convenience wrappers
        function obj = setBatch(obj, id)
            obj = obj.setAnnotation("batch", id);
        end
    
        function obj = setSample(obj, id)
            obj = obj.setAnnotation("sample", id);
        end
    
        function obj = setCluster(obj, labels)
            obj = obj.setAnnotation("cluster", labels);
        end        
    end
    
    methods (Static)
        %% Build from raw matrix
        function obj = fromMatrix(X, g, cellIDs)
            obj = SingleCellExperiment2(X, g, cellIDs);
        end

        %% Convert from SingleCellExperiment
        function obj = fromSCE(sce)
            obj = sce.toSCE2();
        end
        
        %% Load demo dataset
        function obj = demoData(name)
            arguments
                name (1,1) string {mustBeMember(name,["pbmc3k","toy","empty"])}
            end
            
            switch name
                case "pbmc3k"
                    data = load("pbmc3k_demo.mat");
                    if isfield(data, "sce2")
                        obj = data.sce2;
                    elseif isfield(data, "X") && isfield(data, "g")
                        obj = SingleCellExperiment2.fromMatrix(data.X, data.g);
                    else
                        error("pbmc3k_demo.mat must contain sce2 or (X,g).");
                    end
                    
                case "toy"
                    X = poissrnd(1, 50, 10);
                    g = "Gene" + (1:50)';
                    obj = SingleCellExperiment2.fromMatrix(X, g);
                    
                case "empty"
                    obj = SingleCellExperiment2();
            end
            
            obj.metadata.demo = true;
        end
    end
end
