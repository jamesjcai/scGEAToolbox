function [X, g, b, batchid, celltype, filenm] = sc_readh5adfile(filenm)
% Read H5AD file
% https://anndata.readthedocs.io/en/latest/fileformat-prose.html
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
%
% Handles /X stored either as a sparse CSR group (data/indices/indptr, the
% common layout) OR as a dense 2-D dataset (as written by AnnData when X is a
% plain ndarray). Gene/barcode names are read from whatever dataset the AnnData
% "_index" attribute points to, so non-standard index names (e.g. /var/gene,
% /obs/cell instead of /var/_index) are handled as well.

X = [];
g = [];
b = [];
batchid = []; celltype = [];
if nargin < 1 || isempty(filenm)
    [filenm, pathname] = uigetfile( ...
        {'*.h5ad', 'H5AD Files (*.h5ad)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a H5AD file');
    if isequal(filenm, 0), return; end
    filenm = fullfile(pathname, filenm);
end
if exist(filenm, 'file') ~= 2, error('File Not Found.'); end

hinfo = h5info(filenm);
groupNames = strtrim(string(char(hinfo.Groups.Name)));
XisSparse = any(groupNames == "/X");

if XisSparse
    % ---- /X is a sparse CSR group: data / indices / indptr ----------------
    idx = find(groupNames == "/X");

    data = pkg.e_guessh5field(filenm, {'/X/'}, {'data'}, true);
    shapeGroupIdx = idx;  % default: read shape from /X
    rawIdx = [];
    rawXIdx = [];
    if isequal(data(1:5), round(data(1:5)))
        indices = pkg.e_guessh5field(filenm, {'/X/'}, {'indices'}, true);
        indptr = pkg.e_guessh5field(filenm, {'/X/'}, {'indptr'}, true);
    else
        warning('sc_readh5adfile:NormalizedX', ...
            '/X appears transformed/normalized. Attempting to read /raw/X instead.');
        try
            data = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'data'}, true);
            indices = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indices'}, true);
            indptr = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indptr'}, true);
            % Update shape group to /raw/X if it exists
            rawIdx = find(strcmp(strtrim(string(char(hinfo.Groups.Name))), "/raw"));
            if ~isempty(rawIdx)
                rawSubNames = strtrim(string(char(hinfo.Groups(rawIdx).Groups.Name)));
                rawXIdx = find(strcmp(rawSubNames, "/raw/X"));
                if ~isempty(rawXIdx) && ~isempty(hinfo.Groups(rawIdx).Groups(rawXIdx).Attributes)
                    shapeGroupIdx = [];  % signal to use raw group below
                end
            end
        catch ME
            warning('sc_readh5adfile:RawXFailed', ...
                '/raw/X could not be read (%s). Using normalized /X instead.', ME.message);
            indices = pkg.e_guessh5field(filenm, {'/X/'}, {'indices'}, true);
            indptr = pkg.e_guessh5field(filenm, {'/X/'}, {'indptr'}, true);
        end
    end

    if ~isempty(shapeGroupIdx)
        grpAttrs = hinfo.Groups(shapeGroupIdx).Attributes;
    else
        grpAttrs = hinfo.Groups(rawIdx).Groups(rawXIdx).Attributes;
    end
    idx2 = find(strcmp(strtrim(string(char(grpAttrs.Name))), "shape"));
    if isempty(idx2)
        idx2 = find(strcmp(strtrim(string(char(grpAttrs.Name))), "h5sparse_shape"));
    end
    shape = double(grpAttrs(idx2).Value);

    % reconstruct genes-by-cells sparse matrix from CSR (rows = cells)
    if ~isMATLABReleaseOlderThan('R2025a')
        X = spalloc(shape(2), shape(1), length(data), 'single');
    else
        X = spalloc(shape(2), shape(1), length(data));
    end
    for k = 1:length(indptr) - 1
        ix = indptr(k) + 1:indptr(k+1);
        y = indices(ix) + 1;
        X(y, k) = data(ix);
    end
else
    % ---- /X is a dense 2-D dataset ---------------------------------------
    % MATLAB reverses HDF5 dimension order, so h5read of an (n_obs x n_vars)
    % AnnData array returns it as (n_vars x n_obs) = genes-by-cells directly.
    Xd = h5read(filenm, '/X');
    probe = Xd(1:min(5, numel(Xd)));
    if ~isempty(probe) && ~isequal(probe(:), round(probe(:))) && any(groupNames == "/raw")
        warning('sc_readh5adfile:NormalizedX', ...
            '/X appears transformed/normalized. Attempting to read /raw/X instead.');
        try
            Xd = h5read(filenm, '/raw/X');   % dense raw counts, if present
        catch
            % /raw/X absent or sparse; keep the (normalized) dense /X as-is
        end
    end
    X = sparse(double(Xd));
end

% ---- gene names: prefer the dataset named by the var "_index" attribute ---
g = readDataFrameIndex(filenm, '/var');
if isempty(g)
    g = pkg.e_guessh5field(filenm, {'/var/'}, {'_index', 'gene_ids', ...
        'gene_name', 'symbol'}, false);
end
if isempty(g) || isscalar(unique(strlength(g))) % suggesting ENSEMBLE ID
    % Gene IDs look uniform-length (e.g. ENSEMBL); try feature_name for symbols
    gx = pkg.e_guessh5field(filenm, {'/raw/var/feature_name/', ...
        '/var/feature_name/'}, {'categories'}, false);
    if ~isempty(gx)
        g = gx;
    end
end
if isempty(g), warning('sc_readh5adfile:NoGenenames', 'Genename is not assigned.'); end

% ---- barcodes: prefer the dataset named by the obs "_index" attribute -----
b = readDataFrameIndex(filenm, '/obs');
if isempty(b)
    b = pkg.e_guessh5field(filenm, {'/obs/'}, {'_index', 'barcodes', ...
        'cell_id', 'CellID'});
end
if isempty(b), warning('Barcode is not assigned.'); end

try
    batchid = readObsColumn(filenm, 'BatchID');
catch
    % BatchID is optional metadata; absence is not an error
end

try
    celltype = readObsColumn(filenm, 'CellType');
catch
    % CellType is optional metadata; absence is not an error
end

g = deblank(string(g));

end


function names = readDataFrameIndex(h5file, groupPath)
% READDATAFRAMEINDEX Read an AnnData obs/var index regardless of its name.
%   AnnData stores the DataFrame index in a dataset whose name is given by the
%   group's "_index" attribute (commonly "_index", but may be e.g. "gene" or
%   "cell"). This resolves that attribute and reads the corresponding dataset.

    names = string.empty;
    try
        idxName = h5readatt(h5file, groupPath, '_index');
    catch
        idxName = '_index';
    end
    idxName = char(string(idxName));
    try
        raw = h5read(h5file, [groupPath '/' idxName]);
        if ischar(raw) || iscellstr(raw)
            raw = string(raw);
        end
        names = deblank(string(raw));
        names = names(:);
    catch
        names = string.empty;
    end
end


function values = readObsColumn(h5file, colname)
% READOBSCOLUMN Extract a column from AnnData .h5ad obs table.
%   values = READOBSCOLUMN(h5file, colname)
%   h5file : path to .h5ad file
%   colname: string, name of the obs column (e.g. 'celltype')
%
%   Returns either a string array or categorical array.

    obsPath = ['/obs/' colname];

    try
        % Case 1: column stored directly (string array, numeric, etc.)
        values = h5read(h5file, obsPath);
        % Convert to MATLAB string if it's char data
        if ischar(values)
            values = string(values);
        elseif iscellstr(values)
            values = string(values);
        elseif isstring(values)
            % good -- do nothing
        end
        return
    catch
        % If direct read fails, probably categorical
    end

    % Case 2: categorical (codes + categories)
    codesPath = [obsPath '/codes'];
    catsPath  = [obsPath '/categories'];

    try
        codes = h5read(h5file, codesPath);
        cats  = h5read(h5file, catsPath);

        % Convert categories to string
        if ischar(cats) || iscellstr(cats)
            cats = string(cats);
        end

        % AnnData categorical codes are 0-based, MATLAB is 1-based
        values = categorical(codes + 1, 1:numel(cats), cats);
        return
    catch
        values = string.empty;  % column not found in either format
    end
end
