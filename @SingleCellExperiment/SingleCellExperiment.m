classdef SingleCellExperiment < handle & matlab.mixin.Copyable
    properties
        X {mustBeNumeric, mustBeFinite, mustBeNonNan} % counts
        g string % genelist
        s double{mustBeNumeric, mustBeFinite} % cell embeddings
        c % current/active group/class id
        c_cell_cycle_tx % cell cycle string array
        c_cell_type_tx % cell type string array
        c_cluster_id % clustering result vector
        c_batch_id % batch id vector
        c_cell_id % barcode vector
        list_cell_attributes cell % e.g., attributes = {'size',[4,6,2]};
        list_gene_attributes cell % e.g., attributes = {'size',[4,6,2]};
        metadata string
        struct_cell_embeddings = pkg.e_makeembedstruct;
        struct_cell_clusterings = pkg.e_makecluststruct;
        % struct_embed struct
        % struct_clust struct
        % struct_embed struct
        % struct_clust struct
        % c_batch_id string = string.empty(0,1)  % Column vector of strings
        % c_cell_id string = string.empty(0,1)   % Column vector of strings
        % c_cluster_id double = double.empty(0,1) % Column vector of doubles
        % c_batch_id string {mustBeNonmissing}
        % c_cell_id string {mustBeNonmissing}
        % metadata string = string.empty        
    end

    properties (Dependent)
        NumCells
        NumGenes
    end

methods

    function obj = SingleCellExperiment(X, g, s, c)
        if nargin < 1, X = []; end
        % if nargin < 2 || isempty(g), g = pkg.i_num2strcell(size(X, 1), "g"); end
        if nargin < 2 || isempty(g), g = "g"+string((1:size(X, 1))'); end
        if nargin < 3 || isempty(s), s = randn(size(X, 2), 3); end
        if nargin < 4 || isempty(c), c = ones(size(X, 2), 1); end
        assert(size(X, 2) == size(s, 1))
        if ~(size(s, 2) > 2)
            s = [s, zeros(size(X, 2), 1)];
        end
        if ~issparse(X)
            X = sparse(X);
        end

        % Convert to single sparse if supported (R2025a+)
        if isMATLABReleaseOlderThan('R2025a')
            % warning("Single-precision sparse not supported. Keeping double precision.");
        else
            X = single(X); % Works in R2025a+
        end
            
        obj.X = X;
        obj.g = g;
        obj.s = s;
        obj.c = c;
        obj.c_cell_id = string(transpose(1:size(X, 2)));
        obj.c_batch_id = string(ones(size(X, 2), 1));
        obj.c_cluster_id = ones(size(X, 2), 1);
        obj.c_cell_cycle_tx = repmat("undetermined", size(X, 2), 1);
        obj.c_cell_type_tx = repmat("undetermined", size(X, 2), 1);
        obj.metadata = string(sprintf('Created: %s', datetime()));
        % obj.struct_cell_embeddings=struct('tsne',[],'umap',[],'phate',[]);
        %s_embed(1) = struct('id','tsne2d','name','tSNE (2D)','value',[]);
        %s_embed(2) = struct('id','tsne3d','name','tSNE (3D)','value',[]);
        %s_embed(3) = struct('id','umap2d','name','UMAP (2D)','value',[]);
        %s_embed(4) = struct('id','umap3d','name','UMAP (3D)','value',[]);
        %s_embed(5) = struct('id','phate2d','name','PHATE (2D)','value',[]);
        %s_embed(6) = struct('id','phate3d','name','PHATE (3D)','value',[]);
        %s_clust(1) = struct('id','kmeans','name','k-Means','value',[]);
        %s_clust(2) = struct('id','snndpc','name','SNNDPC','value',[]);
        %s_clust(3) = struct('id','sc3','name','SC3','value',[]);
        %obj.struct_embed = s_embed;
        %obj.struct_clust = s_clust;
    end

    function set.s(obj, value)
        % Ensure s is numeric (or you can relax this if needed)
        % validateattributes(value, {'numeric'}, {}, mfilename, 's');
        if size(value, 1) ~= numcells(obj)
            error(['Number of rows in s (%d) must equal ', ...
                   'number of columns in X (%d).'], ...
                   size(value, 1), numcells(obj));
        end
        obj.s = value;
    end

    % function obj = set.g(obj, value)
    %     arguments
    %         obj
    %         value string
    %     end
    %     if numel(value) ~= size(obj.X, 1)
    %         error(['Length of g (%d) must equal ', ...
    %                'number of rows in X (%d).'], ...
    %                numel(value), size(obj.X, 1));
    %     end
    %     obj.g = value;
    % end


    function set.g(obj, value)
        arguments
            obj SingleCellExperiment
            value string {mustBeNonmissing}
        end
        expectedLength = numgenes(obj);
        actualLength = numel(value);
        if actualLength ~= expectedLength
            error('SingleCellExperiment:InvalidGeneList', ...
                'Gene list length (%d) must match number of genes in expression matrix (%d)', ...
                actualLength, expectedLength);
        end
        obj.g = value;
    end    

    function m = get.NumCells(obj)
        m = size(obj.X, 2);
    end

    function set.NumCells(obj, ~)
        fprintf('%s%d\n', 'NumCells is: ', obj.NumCells)
        error('You cannot set NumCells property');
    end

    function m = get.NumGenes(obj)
        m = size(obj.X, 1);
    end

    function set.NumGenes(obj, ~)
        fprintf('%s%d\n', 'NumGenes is: ', obj.NumGenes)
        error('You cannot set NumGenes property');
    end

    function r = numcells(obj)
        r = size(obj.X, 2);
    end

    function r = numgenes(obj)
        r = size(obj.X, 1);
    end

    obj = assigncelltype(obj, speciesid, keepclusterid)
    obj = clustercells(obj, k, methodid, forced, sx)
    obj = embedcells(obj, methodid, forced, usehvgs, ndim, numhvg, whitelist, showwaitbar)
    obj = estimatecellcycle(obj, forced, methodid)
    obj = estimatepotency(obj, speciesid, forced)
    obj = onestepraw2anno(obj, speciesid);
    obj = qcfilterwhitelist(obj, libszcutoff, mtratio, min_cells_nonzero, gnnumcutoff, whitelist)
    obj = sortcells(obj, idx);
    exportToJsonl(obj, outfile, dataset_id, titleText)
   
    function obj = removecells(obj, idx)
        try
            obj.X(:, idx) = [];
        catch ME
            if issparse(obj.X)
                rethrow(ME);
            else
                obj.X = sparse(obj.X);
                obj.X(:, idx) = [];
            end
        end

        obj.s(idx, :) = [];
        obj.c(idx) = [];
        if ~isempty(obj.c_cell_cycle_tx)
            obj.c_cell_cycle_tx(idx) = [];
        end
        if ~isempty(obj.c_cell_type_tx)
            obj.c_cell_type_tx(idx) = [];
        end
        if ~isempty(obj.c_cluster_id)
            obj.c_cluster_id(idx) = [];
        end
        if ~isempty(obj.c_batch_id)
            obj.c_batch_id(idx) = [];
        end
        if ~isempty(obj.c_cell_id)
            obj.c_cell_id(idx) = [];
        end
        for k = 2:2:length(obj.list_cell_attributes)
            obj.list_cell_attributes{k}(idx) = [];
        end

        a = fieldnames(obj.struct_cell_embeddings);
        for k = 1:length(a)
            if ~isempty(obj.struct_cell_embeddings.(a{k}))
                obj.struct_cell_embeddings.(a{k})(idx, :) = [];
            end
        end

        a = fieldnames(obj.struct_cell_clusterings);
        for k = 1:length(a)
            if ~isempty(obj.struct_cell_clusterings.(a{k}))
                obj.struct_cell_clusterings.(a{k})(idx) = [];
            end
        end
        % obj.NumCells=size(obj.X,2);
    end

    function obj = selectcells(obj, idx)
        if islogical(idx)
            if length(idx) == obj.NumCells
                obj = removecells(obj, ~idx);
            else
                error('length(idx)~=sce.NumCells');
            end
        else
            try
                obj.X = obj.X(:, idx);
            catch ME
                if issparse(obj.X)
                    rethrow(ME);
                else
                    obj.X = sparse(obj.X);
                    obj.X = obj.X(:, idx);
                end
            end
            obj.s = obj.s(idx, :);
            obj.c = obj.c(idx);
            if ~isempty(obj.c_cell_cycle_tx)
                obj.c_cell_cycle_tx = obj.c_cell_cycle_tx(idx);
            end
            if ~isempty(obj.c_cell_type_tx)
                obj.c_cell_type_tx = obj.c_cell_type_tx(idx);
            end
            if ~isempty(obj.c_cluster_id)
                obj.c_cluster_id = obj.c_cluster_id(idx);
            end
            if ~isempty(obj.c_batch_id)
                obj.c_batch_id = obj.c_batch_id(idx);
            end
            if ~isempty(obj.c_cell_id)
                obj.c_cell_id = obj.c_cell_id(idx);
            end
            for k = 2:2:length(obj.list_cell_attributes)
                obj.list_cell_attributes{k} = obj.list_cell_attributes{k}(idx);
            end

            a = fieldnames(obj.struct_cell_embeddings);
            for k = 1:length(a)
                if ~isempty(obj.struct_cell_embeddings.(a{k}))
                    obj.struct_cell_embeddings.(a{k}) = ...
                        obj.struct_cell_embeddings.(a{k})(idx, :);
                end
            end

            a = fieldnames(obj.struct_cell_clusterings);
            for k = 1:length(a)
                if ~isempty(obj.struct_cell_clusterings.(a{k}))
                    obj.struct_cell_clusterings.(a{k}) = ...
                        obj.struct_cell_clusterings.(a{k})(idx);
                end
            end
            % obj.NumCells=size(obj.X,2);
        end
    end

    function set.c(obj, tmpc)
        if length(tmpc) ~= numcells(obj)
            error('length(c)~=numcells(sce)');
        else
            obj.c = tmpc;
        end
    end

    function r = title(obj)
        %        r=sprintf('%d x %d\n[genes x cells]',...
        %            size(obj.X,1),size(obj.X,2));
        r = sprintf('%d x %d', ...
            size(obj.X, 1), size(obj.X, 2));
    end

    function obj = qcfilter(obj, libsize, mtratio, min_cells_nonzero)
        % arguments
        %     obj SingleCellExperiment
        %     options.libsize (1,1) double {mustBePositive} = 1000
        %     options.mtratio (1,1) double {mustBeInRange(options.mtratio,0,1)} = 0.15
        %     options.min_cells_nonzero (1,1) double {mustBeNonnegative} = 15
        % end        
        if nargin < 4 || isempty(min_cells_nonzero), min_cells_nonzero = 15; end
        if nargin < 3 || isempty(mtratio), mtratio = 0.15; end
        if nargin < 2 || isempty(libsize), libsize = 1000; end
        %        case 'Relaxed (keep more cells/genes)'
        %            definput = {'500','0.20','10'};
        %        case 'Strigent (keep less cells/genes)'
        %            definput = {'1000','0.15','15'};
        [obj.X, obj.g] = sc_rmdugenes(obj.X, obj.g);
        [~, keptg, keptidxv] = sc_qcfilter(obj.X, obj.g, ...
            libsize, mtratio, ...
            min_cells_nonzero);
        for k = 1:length(keptidxv)
            obj = selectcells(obj, keptidxv{k}); % OK
        end
        [y] = ismember(obj.g, keptg);
        obj.X = obj.X(y, :);
        obj.g = obj.g(y);
    end

    function obj = selectgenes(obj, min_cellnum, nonzero_cutoff)
        if nargin < 3, nonzero_cutoff = 1; end
        if nargin < 2, min_cellnum = 0.01; end
        [tmpX, tmpg, idx] = sc_selectg(obj.X, obj.g, ...
            min_cellnum, nonzero_cutoff);
        obj.X = tmpX;
        obj.g = tmpg;
        try
        for k = 2:2:length(obj.list_gene_attributes)
            obj.list_gene_attributes{k} = obj.list_gene_attributes{k}(idx);
        end
        catch ME
            warning(ME.message);
        end
    end

    function obj = selectkeepgenes(obj, min_countnum, min_cellnum)
        if nargin < 2, min_countnum = 1; end
        if nargin < 3, min_cellnum = 0.01; end
        nc = sum(obj.X >= min_countnum, 2);
        if min_cellnum < 1
            idxkeep1 = nc >= min_cellnum * size(obj.X, 2);
        else
            idxkeep1 = nc >= min_cellnum;
        end
        idxkeep2 = true(size(idxkeep1));
        k = sum(~idxkeep1);
        [~, idx2] = mink(mean(obj.X, 2), k);
        idxkeep2(idx2) = false;
        idxkeep = idxkeep1 | idxkeep2;
        obj.X = obj.X(idxkeep, :);
        obj.g = obj.g(idxkeep);
        try
            for k = 2:2:length(obj.list_gene_attributes)
                obj.list_gene_attributes{k} = obj.list_gene_attributes{k}(idxkeep);
            end
        catch ME
            warning(ME.message);
        end
    end

    function newobj = subsetcopy(obj, idx)
        %SUBSETCOPY Return a copy of the object with only selected cells
        %
        % Usage:
        %   sce1 = sce.subsetcopy(j1);
        %
        newobj = copy(obj);           % make a deep/shallow copy (depending on your class)
        newobj.selectcells(idx);      % run selection on the copy
    end

    function obj = rmmtgenes(obj)
        [tmpX, tmpg, idx] = sc_rmmtgenes(obj.X, obj.g, 'mt-', true);
        if sum(idx) > 0
            obj.X = tmpX;
            obj.g = tmpg;
            try
                for k = 2:2:length(obj.list_gene_attributes)
                    obj.list_gene_attributes{k}(idx) = [];
                end
            catch ME
                warning(ME.message);
            end
        end
    end

    function obj = rmlncrnagenes(obj)
        glist = pkg.i_get_lncrnagenes;
        [idx] = ismember(upper(obj.g), glist);
        obj.X = obj.X(~idx, :);
        obj.g = obj.g(~idx);
        fprintf('%d lncRNA genes found and removed.\n', ...
            sum(idx));
            try
                for k = 2:2:length(obj.list_gene_attributes)
                    obj.list_gene_attributes{k}(idx) = [];
                end
            catch ME
                warning(ME.message);
            end
    end

    function obj = rmribosomalgenes(obj)
        ribog = pkg.i_get_ribosomalgenes;
        [idx] = ismember(upper(obj.g), ribog);
        obj.X = obj.X(~idx, :);
        obj.g = obj.g(~idx);
        fprintf('%d ribosomal genes found and removed.\n', ...
            sum(idx));
            try
                for k = 2:2:length(obj.list_gene_attributes)
                    obj.list_gene_attributes{k}(idx) = [];
                end
            catch ME
                warning(ME.message);
            end
    end

    function obj = rmhemoglobingenes(obj)
        hemog = pkg.i_get_hemoglobingenes;        
        [idx] = ismember(obj.g, hemog);
        obj.X = obj.X(~idx, :);
        obj.g = obj.g(~idx);
        fprintf('%d hemoglobin genes found and removed.\n', ...
            sum(idx));
            try
                for k = 2:2:length(obj.list_gene_attributes)
                    obj.list_gene_attributes{k}(idx) = [];
                end
            catch ME
                warning(ME.message);
            end
    end

    function obj = appendmetainfo(obj, infostr)
        if ~isstring(infostr)
            infostr = string(infostr);
        end
        obj.metadata = [obj.metadata; infostr];
    end

    function obj = setbatchid(obj, id)
        if ischar(id), id = string(id); end
        if isscalar(id)
            obj.c_batch_id = repmat(id, obj.NumCells, 1);
        elseif numel(id) == obj.NumCells
            obj.c_batch_id = id;
        end
    end

    function c_check(obj)
        assert(~isempty(obj.c), 'SCE.C must be defined!');
    end
end
% https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html



end
