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
        % Additional assays measured on THESE cells: CITE-seq ADT, scATAC
        % peaks, and so on. One field per modality; see setModality for the
        % shape of each. A plain struct() default rather than a pkg.e_make*
        % factory, because the field names are open-ended and a literal
        % removes a class-load-time dependency.
        struct_modalities = struct();
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

    obj = assigncelltype(obj, speciesid, keepclusterid, keepold)
    obj = clustercells(obj, k, methodid, forced, sx)
    obj = embedcells(obj, methodid, forced, usehvgs, ndim, numhvg, whitelist, showwaitbar)
    obj = estimatecellcycle(obj, forced, methodid)
    obj = estimatepotency(obj, speciesid, forced)
    obj = onestepraw2anno(obj, speciesid);
    obj = qcfilterwhitelist(obj, libszcutoff, mtratio, min_cells_nonzero, gnnumcutoff, whitelist)
    obj = sortcells(obj, idx);
    exportToJsonl(obj, outfile, dataset_id, titleText)
    sce2 = toSCE2(obj)

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
        obj = i_applyCellIndex(obj, idx, 'delete');
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
            obj = i_applyCellIndex(obj, idx, 'select');
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

    function tf = hasCellAttribute(obj, name)
        % Returns true if a cell attribute with the given name exists.
        tf = any(strcmp(name, obj.list_cell_attributes(1:2:end)));
    end

    function val = getCellAttribute(obj, name)
        % Returns the value vector for a named cell attribute, or [] if absent.
        idx = find(strcmp(name, obj.list_cell_attributes(1:2:end)), 1);
        if isempty(idx)
            val = [];
        else
            val = obj.list_cell_attributes{2*idx};
        end
    end

    function setCellAttribute(obj, name, value)
        % Add or update a named cell attribute.
        % value must be a vector of length NumCells.
        idx = find(strcmp(name, obj.list_cell_attributes(1:2:end)), 1);
        if isempty(idx)
            obj.list_cell_attributes = [obj.list_cell_attributes, {char(name), value}];
        else
            obj.list_cell_attributes{2*idx} = value;
        end
    end

    function removeCellAttribute(obj, name)
        % Remove a named cell attribute if it exists.
        idx = find(strcmp(name, obj.list_cell_attributes(1:2:end)), 1);
        if ~isempty(idx)
            obj.list_cell_attributes([2*idx-1, 2*idx]) = [];
        end
    end

    function tf = hasGeneAttribute(obj, name)
        % Returns true if a gene attribute with the given name exists.
        tf = any(strcmp(name, obj.list_gene_attributes(1:2:end)));
    end

    function val = getGeneAttribute(obj, name)
        % Returns the value vector for a named gene attribute, or [] if absent.
        idx = find(strcmp(name, obj.list_gene_attributes(1:2:end)), 1);
        if isempty(idx)
            val = [];
        else
            val = obj.list_gene_attributes{2*idx};
        end
    end

    function setGeneAttribute(obj, name, value)
        % Add or update a named gene attribute.
        % value must be a vector of length NumGenes.
        idx = find(strcmp(name, obj.list_gene_attributes(1:2:end)), 1);
        if isempty(idx)
            obj.list_gene_attributes = [obj.list_gene_attributes, {char(name), value}];
        else
            obj.list_gene_attributes{2*idx} = value;
        end
    end

    function removeGeneAttribute(obj, name)
        % Remove a named gene attribute if it exists.
        idx = find(strcmp(name, obj.list_gene_attributes(1:2:end)), 1);
        if ~isempty(idx)
            obj.list_gene_attributes([2*idx-1, 2*idx]) = [];
        end
    end

    function obj = setModality(obj, name, X, f, type, present)
        % SETMODALITY Store an assay measured on this object's cells.
        %
        %   sce.setModality("atac", Xpeaks, peakNames)
        %   sce.setModality("prot", Xadt, abNames, "ADT", presentMask)
        %
        %   X must be features-by-cells with EXACTLY NumCells columns, in the
        %   same cell order as sce.X. That is the whole contract: every
        %   cell-axis operation subsets modality columns positionally, so a
        %   width mismatch is silent corruption rather than an error later.
        %   It is rejected here, at the only point where it is still cheap to
        %   diagnose. Use PRESENT to mark cells the assay did not measure -
        %   a cell absent from an ATAC run is not a cell with zero
        %   accessibility, and the two must not be conflated.
        %
        %   X is stored as single, sparse or dense by measured density. Note
        %   that the {mustBeFinite, mustBeNonNan} validators on the X PROPERTY
        %   do not reach a modality, which lives inside a struct: NaN and Inf
        %   are accepted here. That is deliberate - a normalized modality can
        %   legitimately contain them - but it does mean a modality is not
        %   guaranteed finite the way sce.X is.
        name = string(name);
        fieldname = matlab.lang.makeValidName(name);
        if fieldname ~= name
            warning('SingleCellExperiment:ModalityRenamed', ...
                'Modality "%s" stored as "%s" to be a valid field name.', ...
                name, fieldname);
        end

        % Two different names can mangle to the same field ("my mod" and
        % "myMod" both give myMod), which would overwrite one modality with
        % another and leave no trace. Refuse; re-setting the SAME name is
        % still allowed and is the documented way to replace a modality.
        existing = [];
        if isfield(obj.struct_modalities, fieldname)
            existing = obj.struct_modalities.(fieldname);
            if isfield(existing, 'name') && existing.name ~= name
                error('SingleCellExperiment:ModalityNameCollision', ...
                    ['Modality "%s" and the existing "%s" both map to the ', ...
                     'field "%s". Rename one of them.'], ...
                    name, existing.name, fieldname);
            end
        end

        % Replacing a modality keeps whatever the caller did not restate:
        % re-setting the matrix should not silently discard a PRESENT mask
        % or the feature annotations that describe the same features.
        if nargin < 5 || isempty(type)
            if ~isempty(existing), type = existing.type; else, type = "unknown"; end
        end
        if nargin < 6 || isempty(present)
            if ~isempty(existing)
                present = existing.present;
            else
                present = true(1, numcells(obj));
            end
        end

        if ~ismatrix(X)
            error('SingleCellExperiment:ModalityWidth', ...
                ['Modality "%s" must be a 2-D features-by-cells matrix; ', ...
                 'got a %d-D array.'], name, ndims(X));
        end
        if size(X, 2) ~= numcells(obj)
            error('SingleCellExperiment:ModalityWidth', ...
                ['Modality "%s" has %d columns but the object has %d cells. ', ...
                 'A modality must be features-by-cells over the same cells, ', ...
                 'in the same order.'], name, size(X, 2), numcells(obj));
        end
        % 'ABC'(:) is a 3-by-1 char array and string() would make three
        % one-letter feature names, so a single char name must be converted
        % before flattening.
        if ischar(f), f = string(cellstr(f)); end
        f = string(f(:));
        if numel(f) ~= size(X, 1)
            error('SingleCellExperiment:ModalityFeatures', ...
                'Modality "%s" has %d rows but %d feature name(s).', ...
                name, size(X, 1), numel(f));
        end
        present = logical(present(:))';
        if numel(present) ~= numcells(obj)
            error('SingleCellExperiment:ModalityPresent', ...
                'PRESENT has %d elements but the object has %d cells.', ...
                numel(present), numcells(obj));
        end

        m = struct();
        m.name = name;
        m.X = i_storagefor(obj, X);
        m.f = f;
        m.type = string(type);
        m.present = present;
        % Feature annotations describe the features, so they survive a
        % re-set only while the feature axis is unchanged.
        if ~isempty(existing) && isequal(existing.f, f)
            m.featureAnn = existing.featureAnn;
        else
            m.featureAnn = table('Size', [numel(f), 0], 'VariableTypes', {});
        end
        obj.struct_modalities.(fieldname) = m;
    end

    function m = getModality(obj, name)
        % GETMODALITY The named modality struct, or [] if there is none.
        fieldname = matlab.lang.makeValidName(string(name));
        if isfield(obj.struct_modalities, fieldname)
            m = obj.struct_modalities.(fieldname);
        else
            m = [];
        end
    end

    function tf = hasModality(obj, name)
        tf = isfield(obj.struct_modalities, ...
            matlab.lang.makeValidName(string(name)));
    end

    function obj = removeModality(obj, name)
        fieldname = matlab.lang.makeValidName(string(name));
        if isfield(obj.struct_modalities, fieldname)
            obj.struct_modalities = rmfield(obj.struct_modalities, fieldname);
        end
    end

    function names = modalityNames(obj)
        names = string(fieldnames(obj.struct_modalities))';
    end

    function obj = setFeatureAnnotation(obj, name, annName, values)
        % SETFEATUREANNOTATION Annotate a modality's features (peak chr/start,
        % antibody isotype flags, ...). Validated against the feature count on
        % the way in - list_gene_attributes is a flat unvalidated list whose
        % filtering fails silently, and this deliberately does not copy that.
        fieldname = matlab.lang.makeValidName(string(name));
        if ~isfield(obj.struct_modalities, fieldname)
            error('SingleCellExperiment:NoSuchModality', ...
                'No modality named "%s".', name);
        end
        m = obj.struct_modalities.(fieldname);
        values = values(:);
        if numel(values) ~= numel(m.f)
            error('SingleCellExperiment:FeatureAnnLength', ...
                ['Annotation "%s" has %d values but modality "%s" has %d ', ...
                 'features.'], annName, numel(values), name, numel(m.f));
        end
        m.featureAnn.(matlab.lang.makeValidName(string(annName))) = values;
        obj.struct_modalities.(fieldname) = m;
    end

    function v = getFeatureAnnotation(obj, name, annName)
        v = [];
        m = getModality(obj, name);
        if isempty(m), return, end
        annName = matlab.lang.makeValidName(string(annName));
        if ismember(annName, m.featureAnn.Properties.VariableNames)
            v = m.featureAnn.(annName);
        end
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
        obj = i_filterGeneAttributes(obj, y);
    end

    function obj = selectgenes(obj, min_cellnum, nonzero_cutoff)
        if nargin < 3, nonzero_cutoff = 1; end
        if nargin < 2, min_cellnum = 0.01; end
        [tmpX, tmpg, idx] = sc_selectg(obj.X, obj.g, ...
            min_cellnum, nonzero_cutoff);
        obj.X = tmpX;
        obj.g = tmpg;
        obj = i_filterGeneAttributes(obj, idx);
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
        obj = i_filterGeneAttributes(obj, idxkeep);
    end

    function newobj = subsetcopy(obj, idx)
        % SUBSETCOPY Return a copy of the object with only selected cells
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
            obj = i_filterGeneAttributes(obj, ~idx);
        end
    end

    function obj = rmlncrnagenes(obj)
        glist = pkg.i_get_lncrnagenes;
        [idx] = ismember(upper(obj.g), glist);
        obj.X = obj.X(~idx, :);
        obj.g = obj.g(~idx);
        fprintf('%d lncRNA genes found and removed.\n', sum(idx));
        obj = i_filterGeneAttributes(obj, ~idx);
    end

    function obj = rmribosomalgenes(obj)
        ribog = pkg.i_get_ribosomalgenes;
        [idx] = ismember(upper(obj.g), ribog);
        obj.X = obj.X(~idx, :);
        obj.g = obj.g(~idx);
        fprintf('%d ribosomal genes found and removed.\n', sum(idx));
        obj = i_filterGeneAttributes(obj, ~idx);
    end

    function obj = rmhemoglobingenes(obj)
        hemog = pkg.i_get_hemoglobingenes;
        [idx] = ismember(obj.g, hemog);
        obj.X = obj.X(~idx, :);
        obj.g = obj.g(~idx);
        fprintf('%d hemoglobin genes found and removed.\n', sum(idx));
        obj = i_filterGeneAttributes(obj, ~idx);
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

methods (Access = private)
    function obj = i_filterGeneAttributes(obj, idx)
        % Apply logical or subscript index to all gene attributes.
        % idx: logical vector (true = keep) or subscript indices
        try
            for k = 2:2:length(obj.list_gene_attributes)
                obj.list_gene_attributes{k} = obj.list_gene_attributes{k}(idx);
            end
        catch ME
            warning(ME.message);
        end
    end

    function vals = i_getGeneAttributeRows(obj, idx)
        % Capture the values of every gene attribute at the rows in idx,
        % in idx order. Returned as a cell array aligned with the value
        % slots (2:2:end) of list_gene_attributes.
        n = floor(length(obj.list_gene_attributes)/2);
        vals = cell(1, n);
        for k = 1:n
            v = obj.list_gene_attributes{2*k};
            try
                vals{k} = v(idx);
            catch
                vals{k} = [];
            end
        end
    end

    function obj = i_appendGeneAttributeRows(obj, vals)
        % Append previously captured attribute rows (see
        % i_getGeneAttributeRows) to the end of every gene attribute,
        % preserving each attribute's orientation.
        n = floor(length(obj.list_gene_attributes)/2);
        for k = 1:min(n, numel(vals))
            v = obj.list_gene_attributes{2*k};
            a = vals{k};
            if isempty(a), continue, end
            if isrow(v)
                obj.list_gene_attributes{2*k} = [v(:).', a(:).'];
            else
                obj.list_gene_attributes{2*k} = [v(:); a(:)];
            end
        end
    end

    function obj = i_applyCellIndex(obj, idx, mode)
        % Apply index to all cell-level attributes.
        % mode: 'select' keeps idx rows, 'delete' removes idx rows
        if nargin < 3, mode = 'select'; end

        if strcmp(mode, 'select')
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
            obj = i_applyModalityIndex(obj, idx, 'select');
        else  % delete mode
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
            obj = i_applyModalityIndex(obj, idx, 'delete');
        end
    end

    function Y = i_storagefor(~, X)
        % Coerce a modality matrix to the storage the toolbox wants.
        %
        % Single precision, matching what the constructor does to sce.X - but
        % stated as this class's own policy rather than inherited, because
        % there is no set.X and sce.X demonstrably drifts back to double
        % through normalization (sc_norm returns double).
        %
        % Sparsity is chosen by MEASURED density, not by fiat. Sparse single
        % costs 12 bytes per nonzero (4 data + 8 index) against 4 bytes per
        % element dense, so the break-even is one third. ATAC peaks at ~2%
        % belong sparse; a CITE-seq ADT panel is often ~90% detected, where
        % sparse costs about 1.8x more than dense.
        % Single-precision SPARSE only exists from R2025a, which is why the
        % constructor guards its own cast the same way. Dense single is fine
        % on every release, so on older ones a sparse modality stays double
        % rather than failing outright.
        singleSparse = ~isMATLABReleaseOlderThan('R2025a');

        density = nnz(X) / max(numel(X), 1);
        if density > 1/3
            Y = single(full(X));
        elseif singleSparse
            if ~issparse(X), X = sparse(double(X)); end
            Y = single(X);
        else
            Y = sparse(double(X));
        end
    end

    function obj = i_applyModalityIndex(obj, idx, mode)
        % Subset every modality's COLUMNS with the same index the rest of the
        % cell axis just took. Note this covers reordering as well as
        % dropping, because sortcells routes through i_applyCellIndex too.
        %
        % The whole struct is rebuilt in a local and assigned back ONCE.
        % qcfilter replays its accumulated index chain through selectcells
        % several times per call, so a property get-modify-set per modality
        % per replay would copy the matrices far more often than necessary.
        if isempty(fieldnames(obj.struct_modalities)), return, end

        mods = obj.struct_modalities;
        names = fieldnames(mods);
        for k = 1:numel(names)
            m = mods.(names{k});
            if strcmp(mode, 'select')
                m.X = m.X(:, idx);
                m.present = m.present(idx);
            else
                m.X(:, idx) = [];
                m.present(idx) = [];
            end
            mods.(names{k}) = m;
        end
        obj.struct_modalities = mods;
    end
end

% https://www.mathworks.com/help/matlab/matlab_oop/example-representing-structured-data.html

end
