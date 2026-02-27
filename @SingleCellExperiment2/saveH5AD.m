function saveH5AD(obj, filename)
    % Save SingleCellExperiment2 into .h5ad format
    %
    % Features:
    %   - Dense or sparse X
    %   - obs/var with categorical encoding
    %   - obs._index = cell IDs
    %   - var._index = gene IDs
    %   - uns (metadata), obsm (embeddings)

    arguments
        obj
        filename (1,1) string
    end

    if isfile(filename)
        delete(filename);
    end

    %% --- Save X ---
    if issparse(obj.X)
        [i,j,v] = find(obj.X);
        nrow = size(obj.X,1);
        ncol = size(obj.X,2);

        % CSR conversion
        indptr = zeros(nrow+1,1);
        [~,order] = sortrows([i,j]);
        i = i(order); j = j(order); v = v(order);

        for r = 1:nrow
            indptr(r+1) = indptr(r) + nnz(i==r);
        end
        indices = j-1; % 0-based

        h5create(filename,'/X/data',size(v),'Datatype','double');
        h5write(filename,'/X/data',v);

        h5create(filename,'/X/indices',size(indices),'Datatype','int32');
        h5write(filename,'/X/indices',int32(indices));

        h5create(filename,'/X/indptr',size(indptr),'Datatype','int32');
        h5write(filename,'/X/indptr',int32(indptr));

        h5create(filename,'/X/shape',[2 1],'Datatype','int64');
        h5write(filename,'/X/shape',int64([nrow;ncol]));
    else
        h5create(filename,'/X',size(obj.X),'Datatype','double');
        h5write(filename,'/X',obj.X);
    end

    %% --- Save obs ---
    if ~isempty(obj.cellAnn)
        varnames = obj.cellAnn.Properties.VariableNames;

        % Write each column
        for k = 1:numel(varnames)
            col = obj.cellAnn.(varnames{k});
            path = "/obs/" + varnames{k};

            if isstring(col) || iscategorical(col)
                cats = unique(col,'stable');
                codes = zeros(numel(col),1,'int32');
                for ci = 1:numel(cats)
                    codes(col==cats(ci)) = ci-1;
                end
                h5create(filename, path, size(codes),'Datatype','int32');
                h5write(filename, path, codes);

                h5create(filename, path + "/categories", size(cats));
                h5write(filename, path + "/categories", cellstr(cats));
            else
                h5create(filename, path, size(col));
                h5write(filename, path, col);
            end
        end

        % Save _index (cell IDs / rownames)
        if ~isempty(obj.cellAnn.Properties.RowNames)
            idx = obj.cellAnn.Properties.RowNames;
        else
            idx = strcat("cell", string((1:height(obj.cellAnn))'));
        end
        h5create(filename,"/obs/_index",size(idx));
        h5write(filename,"/obs/_index",cellstr(idx));
    end

    %% --- Save var ---
    if ~isempty(obj.geneAnn)
        fnames = fieldnames(obj.geneAnn);
        for k = 1:numel(fnames)
            col = obj.geneAnn.(fnames{k});
            path = "/var/" + fnames{k};

            if isstring(col) || iscategorical(col)
                cats = unique(col,'stable');
                codes = zeros(numel(col),1,'int32');
                for ci = 1:numel(cats)
                    codes(col==cats(ci)) = ci-1;
                end
                h5create(filename, path, size(codes),'Datatype','int32');
                h5write(filename, path, codes);

                h5create(filename, path + "/categories", size(cats));
                h5write(filename, path + "/categories", cellstr(cats));
            else
                h5create(filename, path, size(col));
                h5write(filename, path, col);
            end
        end

        % Save _index (gene IDs)
        if isfield(obj.geneAnn,"names")
            idx = obj.geneAnn.names;
        else
            idx = strcat("gene", string((1:height(struct2table(obj.geneAnn))))');
        end
        h5create(filename,"/var/_index",size(idx));
        h5write(filename,"/var/_index",cellstr(idx));
    end

    %% --- Save uns ---
    if ~isempty(obj.metadata)
        fnames = fieldnames(obj.metadata);
        for k = 1:numel(fnames)
            val = obj.metadata.(fnames{k});
            path = "/uns/" + fnames{k};
            h5create(filename, path, size(val));
            h5write(filename, path, val);
        end
    end

    %% --- Save obsm ---
    if ~isempty(obj.embeddings)
        fnames = fieldnames(obj.embeddings);
        for k = 1:numel(fnames)
            val = obj.embeddings.(fnames{k});
            path = "/obsm/" + fnames{k};
            h5create(filename, path, size(val));
            h5write(filename, path, val);
        end
    end
end