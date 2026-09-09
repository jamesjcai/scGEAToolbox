function [succeeded] = sc_sce2h5ad(sce, filename, verbose, ~)
%SC_SCE2H5AD Write a SingleCellExperiment to an H5AD file.
%
%   succeeded = SC_SCE2H5AD(sce, filename) writes SCE to an AnnData H5AD
%   file that scanpy, anndata and Seurat can open. This is native MATLAB;
%   no Python is involved.
%
%   AnnData is cells-by-genes and SCE is genes-by-cells, so the matrix is
%   transposed on the way out. That costs nothing here: MATLAB stores a
%   sparse matrix column by column, and the column-wise layout of a
%   genes-by-cells matrix is exactly the row-wise (CSR) layout AnnData wants
%   for the cells-by-genes transpose, so the same three arrays serve both.
%
%   What is written:
%     /X              counts, CSR sparse or dense according to density
%     /obs            cell barcodes as the index, plus batch id, cell type,
%                     cluster id, cell cycle phase and every entry of
%                     SCE.LIST_CELL_ATTRIBUTES that has one value per cell
%     /var            gene names as the index, plus SCE.LIST_GENE_ATTRIBUTES
%     /obsm           every embedding in SCE.STRUCT_CELL_EMBEDDINGS, named
%                     X_umap, X_tsne and so on, plus X_sce for SCE.S
%     /uns            the SCE metadata string
%
%   INPUTS:
%     sce      - a SingleCellExperiment.
%     filename - destination path. Prompts when omitted.
%     verbose  - print progress (default true).
%     pe       - former Python environment argument, accepted and ignored so
%                that existing calls keep working. Use RUN.PY_WRITEH5AD
%                directly if the Python writer is still wanted.
%
%   OUTPUT:
%     succeeded - true when the file was written.
%
%   See also SC_READH5ADFILE, SC_SCE2HDF5, RUN.PY_WRITEH5AD.

if nargin < 3 || isempty(verbose), verbose = true; end

succeeded = false;
if nargin < 2 || isempty(filename)
    [filename, pathname] = uiputfile({'*.h5ad'; '*.*'}, 'Save as');
    if ~ischar(filename) && ~isstring(filename), return; end
    filename = fullfile(pathname, filename);
end
filename = char(filename);

if isfile(filename)
    delete(filename);
end

try
    i_writeh5ad(sce, filename, verbose);
    succeeded = true;
catch ME
    if verbose
        fprintf(2, 'sc_sce2h5ad failed: %s\n', ME.message);
    end
    if isfile(filename)
        % A partial file is worse than none: it looks like a result and
        % fails somewhere downstream instead of here.
        delete(filename);
    end
    rethrow(ME);
end

if verbose
    d = dir(filename);
    fprintf('Wrote %s (%d cells x %d genes, %.1f MB)\n', filename, ...
        sce.NumCells, sce.NumGenes, d.bytes/1e6);
end

end


function i_writeh5ad(sce, filename, verbose)

numCells = sce.NumCells;
numGenes = sce.NumGenes;

% ---- /X ---------------------------------------------------------------
X = sce.X;
if ~issparse(X)
    density = nnz(X)/max(numel(X), 1);
    if density < 0.5
        X = sparse(X);
    end
end

if issparse(X)
    % The CSC arrays of the genes-by-cells matrix are the CSR arrays of its
    % cells-by-genes transpose, which is the orientation AnnData stores.
    % SCE.X is often single, and single loses whole integers above 2^24, so
    % a reduction over the values themselves cannot be trusted at scale: on
    % a 7.4M-nonzero matrix MATLAB's own sum(X, "all") was out by 305,700
    % against the exact total. FIND and the count of nonzeros per column are
    % safe -- SUM over a logical returns double -- and the values are cast
    % before they are written.
    [geneIdx, ~, values] = find(X);
    indptr = full(cumsum([0; sum(X ~= 0, 1).']));
    i_writearray(filename, '/X/data', double(values(:)));
    i_writearray(filename, '/X/indices', int32(geneIdx(:) - 1));
    i_writearray(filename, '/X/indptr', int64(indptr(:)));
    i_setattrs(filename, '/X', "csr_matrix", "0.1.0");
    h5writeatt(filename, '/X', 'shape', int64([numCells; numGenes]));
else
    % HDF5 reverses the dimension order of a MATLAB array, so writing the
    % genes-by-cells matrix produces a cells-by-genes dataset on disk.
    h5create(filename, '/X', [numGenes, numCells], 'Datatype', 'double');
    h5write(filename, '/X', full(double(X)));
    i_setattrs(filename, '/X', "array", "0.2.0");
end
if verbose
    if issparse(X)
        fprintf('  /X csr_matrix, %d nonzeros\n', nnz(X));
    else
        fprintf('  /X dense, %d cells x %d genes\n', numCells, numGenes);
    end
end

% ---- /obs -------------------------------------------------------------
barcodes = sce.c_cell_id;
if isempty(barcodes) || numel(barcodes) ~= numCells
    barcodes = "cell" + string((1:numCells).');
end
obs = i_columnstruct();
obs = i_addcolumn(obs, 'batch_id', sce.c_batch_id, numCells);
obs = i_addcolumn(obs, 'cell_type', sce.c_cell_type_tx, numCells);
obs = i_addcolumn(obs, 'cluster_id', sce.c_cluster_id, numCells);
obs = i_addcolumn(obs, 'cell_cycle', sce.c_cell_cycle_tx, numCells);
obs = i_addcolumn(obs, 'active_group', sce.c, numCells);
obs = i_addlistattributes(obs, sce.list_cell_attributes, numCells);
i_writeframe(filename, '/obs', i_uniquestrings(string(barcodes(:))), obs);

% ---- /var -------------------------------------------------------------
genes = string(sce.g(:));
if numel(genes) ~= numGenes
    genes = "gene" + string((1:numGenes).');
end
vars = i_columnstruct();
vars = i_addlistattributes(vars, sce.list_gene_attributes, numGenes);
i_writeframe(filename, '/var', i_uniquestrings(genes), vars);

% ---- /obsm ------------------------------------------------------------
% AnnData draws anything under obsm with the X_ prefix, so name the
% embeddings the way scanpy expects rather than the way SCE stores them.
i_creategroup(filename, '/obsm');
i_setattrs(filename, '/obsm', "dict", "0.1.0");
written = strings(0, 1);
embeddings = sce.struct_cell_embeddings;
if isstruct(embeddings)
    names = string(fieldnames(embeddings));
    for k = 1:numel(names)
        value = embeddings.(names(k));
        if isnumeric(value) && ~isempty(value) && size(value, 1) == numCells
            name = "X_" + lower(names(k));
            i_writematrix(filename, "/obsm/" + name, double(value));
            written(end+1) = name; %#ok<AGROW>
        end
    end
end
if ~isempty(sce.s) && size(sce.s, 1) == numCells && ~any(written == "X_sce")
    i_writematrix(filename, '/obsm/X_sce', double(sce.s));
end

% ---- empty containers AnnData expects ---------------------------------
for group = ["/layers", "/obsp", "/varp", "/varm", "/uns"]
    i_creategroup(filename, char(group));
    i_setattrs(filename, char(group), "dict", "0.1.0");
end
if ~isempty(sce.metadata)
    i_writestringarray(filename, '/uns/metadata', ...
        strjoin(string(sce.metadata(:)).', ' '));
end

% ---- root -------------------------------------------------------------
i_setattrs(filename, '/', "anndata", "0.1.0");

end


function s = i_columnstruct()
s = struct('name', {}, 'value', {});
end


function s = i_addcolumn(s, name, value, expected)
% Keep a column only when it has one entry per cell or gene. SCE fields are
% often empty or stale, and a short column would corrupt the frame.
if isempty(value) || numel(value) ~= expected
    return;
end
if iscellstr(value) || ischar(value) %#ok<ISCLSTR>
    value = string(value);
end
if isstring(value) || iscategorical(value)
    value = string(value(:));
elseif islogical(value)
    value = double(value(:));
elseif isnumeric(value)
    value = double(value(:));
else
    return;
end
s(end+1) = struct('name', name, 'value', value);
end


function s = i_addlistattributes(s, list, expected)
% LIST_CELL_ATTRIBUTES and LIST_GENE_ATTRIBUTES hold {name, value} pairs.
if isempty(list) || ~iscell(list)
    return;
end
for k = 1:2:numel(list) - 1
    name = list{k};
    if ~(ischar(name) || isstring(name))
        continue;
    end
    name = matlab.lang.makeValidName(string(name));
    if any(string({s.name}) == name)
        continue;
    end
    s = i_addcolumn(s, char(name), list{k+1}, expected);
end
end


function i_writeframe(filename, group, index, columns)
% An AnnData dataframe is a group whose attributes name the index dataset
% and fix the column order; the columns themselves are plain datasets.
i_creategroup(filename, group);
i_writestringarray(filename, [group '/_index'], index);
for k = 1:numel(columns)
    path = [group '/' columns(k).name];
    if isstring(columns(k).value)
        i_writestringarray(filename, path, columns(k).value);
    else
        i_writearray(filename, path, columns(k).value);
    end
end
i_setattrs(filename, group, "dataframe", "0.2.0");
i_writestrattr(filename, group, '_index', "_index", true);
if isempty(columns)
    i_writestrattr(filename, group, 'column-order', strings(0, 1), false);
else
    i_writestrattr(filename, group, 'column-order', ...
        string({columns.name}).', false);
end
end


function i_writearray(filename, path, value)
h5create(filename, path, numel(value), 'Datatype', class(value));
h5write(filename, path, value(:));
i_setattrs(filename, path, "array", "0.2.0");
end


function i_writematrix(filename, path, value)
% Written transposed, because HDF5 reverses the dimension order of a MATLAB
% array and obsm entries must come back as nCells-by-nDims.
path = char(path);
h5create(filename, path, [size(value, 2), size(value, 1)], 'Datatype', 'double');
h5write(filename, path, value.');
i_setattrs(filename, path, "array", "0.2.0");
end


function i_writestringarray(filename, path, value)
h5create(filename, path, numel(value), 'Datatype', 'string');
h5write(filename, path, string(value(:)));
i_setattrs(filename, path, "string-array", "0.2.0");
end


function i_setattrs(filename, path, encodingType, encodingVersion)
i_writestrattr(filename, path, 'encoding-type', encodingType, true);
i_writestrattr(filename, path, 'encoding-version', encodingVersion, true);
end


function i_writestrattr(filename, path, name, value, asScalar)
% Write a string attribute through the low-level interface.
%
% H5WRITEATT turns a string scalar into a one-element *array*, and AnnData
% reads encoding-type, encoding-version and _index as scalars: given an
% array it fails in its type registry with "unhashable type: numpy.ndarray",
% nowhere near the attribute that caused it. The same applies at the other
% end, where an empty column-order has to be a zero-length array rather than
% the null dataspace H5WRITEATT produces, or the dataframe reader raises on
% list(). Neither shape is reachable through the high-level interface.
value = string(value);
fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
closeFile = onCleanup(@() H5F.close(fid));
oid = H5O.open(fid, path, 'H5P_DEFAULT');
closeObject = onCleanup(@() H5O.close(oid));

typeId = H5T.copy('H5T_C_S1');
closeType = onCleanup(@() H5T.close(typeId));
H5T.set_size(typeId, 'H5T_VARIABLE');
H5T.set_cset(typeId, H5ML.get_constant_value('H5T_CSET_UTF8'));

if asScalar
    spaceId = H5S.create('H5S_SCALAR');
else
    spaceId = H5S.create_simple(1, numel(value), numel(value));
end
closeSpace = onCleanup(@() H5S.close(spaceId));

% The file is deleted before writing and every attribute is set once, so a
% name collision here would be a bug in this function and should surface.
attrId = H5A.create(oid, name, typeId, spaceId, 'H5P_DEFAULT');
closeAttr = onCleanup(@() H5A.close(attrId));
if asScalar
    H5A.write(attrId, typeId, {char(value)});
elseif ~isempty(value)
    % A zero-length attribute is fully described by its dataspace; asking
    % H5A.write to store no elements fails inside the library.
    H5A.write(attrId, typeId, cellstr(value(:)).');
end
end


function i_creategroup(filename, group)
% There is no h5creategroup, and writing an attribute to a path that does
% not exist fails, so make the group through the low-level interface.
if isfile(filename)
    fid = H5F.open(filename, 'H5F_ACC_RDWR', 'H5P_DEFAULT');
else
    fid = H5F.create(filename, 'H5F_ACC_TRUNC', 'H5P_DEFAULT', 'H5P_DEFAULT');
end
closeFile = onCleanup(@() H5F.close(fid));
plist = 'H5P_DEFAULT';
gid = H5G.create(fid, group, plist, plist, plist);
H5G.close(gid);
end


function s = i_uniquestrings(s)
% AnnData requires unique obs and var names; duplicated gene symbols are
% common enough that silently writing them would break the reader.
s = string(s(:));
s(ismissing(s)) = "";
if numel(unique(s)) ~= numel(s)
    s = string(matlab.lang.makeUniqueStrings(cellstr(s)));
end
end
