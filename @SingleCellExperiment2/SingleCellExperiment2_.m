classdef SingleCellExperiment2_ < handle & matlab.mixin.Copyable
    %SINGLECELLEXPERIMENT Container for single-cell RNA-seq experiments.
    %
    %   Stores gene expression matrix (sparse), cell and gene annotations,
    %   embeddings, clusterings, and metadata.
    %
    %   Example:
    %       sce2 = SingleCellExperiment2(X, g);
    % https://chatgpt.com/share/68afc5eb-6eac-8005-b211-66c0a818ebf9

    properties
        % Core data
        X {mustBeNumeric, mustBeFinite} % sparse counts [genes x cells]

        % Annotations
        cellAnn struct   % cell annotations (cluster, type, cycle, batch, ID, QC)
        geneAnn struct   % gene annotations (names + others)

        % Analysis results
        embeddings struct    % e.g. tsne, umap
        clusterings struct   % e.g. kmeans, sc3

        % Metadata
        metadata struct      % creation info, notes, version, etc.        
    end

    properties (Dependent, SetAccess=private)
        NumCells
        NumGenes
    end

    methods
        %% ---- Constructor ----
        function obj = SingleCellExperiment2_(X, g)
            arguments
                X double {mustBeFinite,mustBeNonNan} = []
                g string = []
            end
            if isempty(X)
                return
            end

            if ~issparse(X), X = sparse(X); end
            if ~isMATLABReleaseOlderThan('R2025a')
                X = single(X);  % single-precision sparse if supported
            end
            obj.X = X;

            if isempty(g)
                g = "g" + string((1:size(X,1))');
            end

            % initialize annotations
            obj.geneAnn = struct(); 
            obj.geneAnn.names = g;

            obj.cellAnn = struct(); 
            obj.cellAnn.ID    = string(1:size(X,2));
            obj.cellAnn.batch = repmat("1", obj.NumCells,1);
            obj.cellAnn.cluster = ones(obj.NumCells,1);
            obj.cellAnn.type  = repmat("undetermined", obj.NumCells,1);
            obj.cellAnn.cycle = repmat("undetermined", obj.NumCells,1);

            % QC defaults
            obj.cellAnn.totalCounts   = full(sum(obj.X,1))';
            obj.cellAnn.detectedGenes = sum(obj.X>0,1)';

            % initialize analysis slots
            obj.embeddings = struct();
            obj.clusterings = struct();

            % metadata
            obj.metadata = struct( ...
                'created', datetime(), ...
                'version', "1.0", ...
                'notes', "" );
        end

        %% ---- Dependent properties ----
        function n = get.NumCells(obj)
            n = size(obj.X,2);
        end
        function n = get.NumGenes(obj)
            n = size(obj.X,1);
        end

        %% ---- Annotation helpers ----
        function obj = setBatch(obj, id)
            if isscalar(id)
                obj.cellAnn.batch = repmat(string(id), obj.NumCells,1);
            elseif numel(id) == obj.NumCells
                obj.cellAnn.batch = string(id);
            else
                error('Batch ID length mismatch');
            end
        end

        %% ---- Selection / subsetting ----
        function obj = subsetCells(obj, idx)
            obj.X = obj.X(:, idx);
            f = fieldnames(obj.cellAnn);
            for k = 1:numel(f)
                if ~isempty(obj.cellAnn.(f{k}))
                    obj.cellAnn.(f{k}) = obj.cellAnn.(f{k})(idx,:);
                end
            end
            fn = fieldnames(obj.embeddings);
            for k = 1:numel(fn)
                if ~isempty(obj.embeddings.(fn{k}))
                    obj.embeddings.(fn{k}) = obj.embeddings.(fn{k})(idx,:);
                end
            end
            fn = fieldnames(obj.clusterings);
            for k = 1:numel(fn)
                if ~isempty(obj.clusterings.(fn{k}))
                    obj.clusterings.(fn{k}) = obj.clusterings.(fn{k})(idx);
                end
            end
        end

        function obj = subsetGenes(obj, idx)
            obj.X = obj.X(idx,:);
            obj.geneAnn.names = obj.geneAnn.names(idx);
            % extend for other gene-level annotations if needed
        end

        %% ---- QC / cleanup methods ----
        function obj = rmEmptyCells(obj)
            keep = sum(obj.X,1) > 0;
            obj.subsetCells(keep);
        end

        function obj = rmEmptyGenes(obj)
            keep = sum(obj.X,2) > 0;
            obj.subsetGenes(keep);
        end

        function obj = rmMitoGenes(obj, mitogenes)
            % Remove mitochondrial genes if supplied
            [~,idx] = intersect(upper(obj.geneAnn.names), upper(mitogenes));
            keep = true(obj.NumGenes,1);
            keep(idx) = false;
            obj.subsetGenes(keep);
        end

        function obj = rmMTGenes(obj)
            % Auto-remove genes starting with "MT-"
            ismt = startsWith(upper(obj.geneAnn.names), "MT-");
            obj.subsetGenes(~ismt);
        end

        function obj = qcFilter(obj, varargin)
            % Flexible QC filter for cells
            p = inputParser;
            addParameter(p,'MinCounts',0,@isscalar);
            addParameter(p,'MaxCounts',inf,@isscalar);
            addParameter(p,'MinGenes',0,@isscalar);
            addParameter(p,'MaxGenes',inf,@isscalar);
            parse(p,varargin{:});
            args = p.Results;

            totalCounts = full(sum(obj.X,1))';
            detectedGenes = sum(obj.X>0,1)';

            keep = (totalCounts >= args.MinCounts) & ...
                   (totalCounts <= args.MaxCounts) & ...
                   (detectedGenes >= args.MinGenes) & ...
                   (detectedGenes <= args.MaxGenes);

            obj.subsetCells(keep);
        end

        %% ---- Metadata helpers ----
        function obj = addNote(obj, str)
            obj.metadata.notes = [obj.metadata.notes; string(str)];
        end

        function obj = addMetadata(obj, note)
            if ~isfield(obj.metadata, "upgradeNotes")
                obj.metadata.upgradeNotes = {};
            end
            obj.metadata.upgradeNotes{end+1} = note;
        end

        function out = g(obj)
            out = obj.geneAnn.names;
        end

        function out = s(obj)
            out = obj.embeddings.('tsne2d');
        end

        %% ---- Display ----
        function disp(obj)
            fprintf('SingleCellExperiment object\n');
            fprintf('  %d genes x %d cells\n', obj.NumGenes, obj.NumCells);
            if isfield(obj.metadata,'created')
                fprintf('  Created: %s\n', string(obj.metadata.created));
            end
            if isfield(obj.metadata,'notes') && ~isempty(obj.metadata.notes)
                fprintf('  Notes: %s\n', strjoin(obj.metadata.notes,"; "));
            end
            fprintf('\nCell annotations: %s\n', strjoin(fieldnames(obj.cellAnn), ", "));
            fprintf('Gene annotations: %s\n', strjoin(fieldnames(obj.geneAnn), ", "));
        end
    end

    methods (Static)
        function obj = fromMatrix(X, g)
            %FROMMATRIX Create a SingleCellExperiment2 from matrix + gene names
            %
            %   obj = SingleCellExperiment2.fromMatrix(X, g)
            %
            %   X : genes x cells expression matrix
            %   g : gene names (string/cell array, length = #genes)
    
            if nargin < 2
                error("fromMatrix requires X (matrix) and g (gene names).");
            end
    
            % Create object
            obj = SingleCellExperiment2_();
    
            % Assign core data
            obj.X = X;
    
            % Ensure geneAnn struct initialized
            obj.geneAnn = struct();
            obj.geneAnn.names = g(:);
    
            % Ensure cellAnn struct initialized
            obj.cellAnn = struct();
            obj.cellAnn.labels = [];
            obj.cellAnn.qualityScores = [];
    
            % Initialize embeddings/clusterings
            obj.embeddings = struct();
            obj.clusterings = struct();
    
            % Metadata
            obj.metadata = struct();
            obj.metadata.source = "fromMatrix";
            obj.metadata.created = datetime;
        end
    
        function obj = demoData(name)
            %DEMODATA Load a demonstration dataset
            %
            %   obj = SingleCellExperiment2.demoData("pbmc3k")
        
            arguments
                name (1,1) string {mustBeMember(name, ["pbmc3k","toy","empty"])}
            end
        
            switch name
                case "pbmc3k"
                    data = load("pbmc3k_demo.mat");
                    
                    if isfield(data, "sce2")
                        % Already a SingleCellExperiment2 object
                        obj = data.sce2;
        
                    elseif isfield(data, "X") && isfield(data, "g")
                        % Build from raw matrix + gene names
                        obj = SingleCellExperiment2_.fromMatrix(data.X, data.g);
        
                    else
                        error("pbmc3k_demo.mat must contain either sce2 or (X,g).");
                    end
        
                case "toy"
                    % Small toy dataset
                    X = poissrnd(1, 50, 10);          % 50 genes × 10 cells
                    g = "Gene" + (1:50)';
                    obj = SingleCellExperiment2_.fromMatrix(X, g);
        
                case "empty"
                    obj = SingleCellExperiment2_();    % returns empty object
                    obj.metadata.source = "emptyDemo";
            end
        
            % Add demo flag
            obj.metadata.demo = true;
        end
    end
    
end