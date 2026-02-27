function obj = loadH5AD(filename)
    % Load a .h5ad file into SingleCellExperiment2
    %
    %   obj = SingleCellExperiment2.loadH5AD("file.h5ad")
    %
    % Supports:
    %   - Dense and CSR sparse /X
    %   - Categorical obs/var decoding
    %   - Metadata (/uns)
    %   - Embeddings (/obsm)

    arguments
        filename (1,1) string
    end

    if ~isfile(filename)
        error("File not found: %s", filename);
    end

    %% --- Load X ---
    X = [];
    try
        X = h5read(filename, '/X');
    catch
        % Try Sparse CSR
        try
            data    = h5read(filename, '/X/data');
            indices = h5read(filename, '/X/indices');
            indptr  = h5read(filename, '/X/indptr');
            shape   = h5read(filename, '/X/shape');

            nrow = shape(1);
            ncol = shape(2);

            rows = zeros(numel(indices),1);
            for i = 1:numel(indptr)-1
                start_idx = indptr(i)+1;
                stop_idx  = indptr(i+1);
                if stop_idx >= start_idx
                    rows(start_idx:stop_idx) = i;
                end
            end
            cols = indices + 1; % Python -> MATLAB
            X = sparse(rows, cols, data, nrow, ncol);

        catch
            error("Could not read /X dataset. Not dense or CSR sparse.");
        end
    end

    %% --- Load obs (cell annotations) ---
    cellAnn = table;
    try
        info = h5info(filename, '/obs');
        for k = 1:numel(info.Datasets)
            dsname = info.Datasets(k).Name;

            % check if categorical
            cats_path = "/obs/" + dsname + "/categories";
            if any(strcmp({info.Datasets.Name}, dsname)) && ...
               h5exists(filename, cats_path)
                % integer codes
                codes = h5read(filename, "/obs/" + dsname) + 1; % 0-based
                cats  = string(h5read(filename, cats_path));
                val   = cats(codes);
            else
                val = h5read(filename, "/obs/" + dsname);
                if isstring(val) || ischar(val)
                    val = string(val);
                end
            end

            if strcmp(dsname, "_index")
                cellAnn.Properties.RowNames = cellstr(val);
            else
                cellAnn.(dsname) = val;
            end
        end
    catch
        % no obs
    end

    %% --- Load var (gene annotations) ---
    geneAnn = table;
    try
        info = h5info(filename, '/var');
        for k = 1:numel(info.Datasets)
            dsname = info.Datasets(k).Name;

            cats_path = "/var/" + dsname + "/categories";
            if h5exists(filename, cats_path)
                codes = h5read(filename, "/var/" + dsname) + 1;
                cats  = string(h5read(filename, cats_path));
                val   = cats(codes);
            else
                val = h5read(filename, "/var/" + dsname);
                if isstring(val) || ischar(val)
                    val = string(val);
                end
            end

            if strcmp(dsname, "_index")
                geneAnn.names = val;
            else
                geneAnn.(dsname) = val;
            end
        end
    catch
        % no var
    end

    %% --- Load uns (metadata) ---
    metadata = struct;
    try
        info = h5info(filename, '/uns');
        for k = 1:numel(info.Datasets)
            dsname = info.Datasets(k).Name;
            metadata.(dsname) = h5read(filename, "/uns/" + dsname);
        end
    catch
        % no uns
    end

    %% --- Load obsm (embeddings) ---
    embeddings = struct;
    try
        info = h5info(filename, '/obsm');
        for k = 1:numel(info.Datasets)
            dsname = info.Datasets(k).Name;
            embeddings.(dsname) = h5read(filename, "/obsm/" + dsname);
        end
    catch
        % no obsm
    end

    %% --- Construct object ---
    obj = SingleCellExperiment2(X, geneAnn, cellAnn);
    obj.metadata   = metadata;
    obj.embeddings = embeddings;
end


%% --- Helper: check existence of dataset ---
function tf = h5exists(filename, path)
    tf = false;
    try
        h5info(filename, path);
        tf = true;
    catch
    end
end
