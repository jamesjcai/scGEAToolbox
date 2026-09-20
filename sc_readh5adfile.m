function [X, g, b, batchid, celltype, filenm] = sc_readh5adfile(filenm)
% Read H5AD file
% https://anndata.readthedocs.io/en/latest/fileformat-prose.html
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices
%
% Handles /X stored either as a sparse group (data/indices/indptr, in CSR or
% CSC layout) OR as a dense 2-D dataset (as written by AnnData when X is a
% plain ndarray). Gene/barcode names are read from whatever dataset the AnnData
% "_index" attribute points to, so non-standard index names (e.g. /var/gene,
% /obs/cell instead of /var/_index) are handled as well. Batch and cell type
% columns are matched against obs ignoring case and separators, so "CellType",
% "cell_type" and "cell type" are all recognized.

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
    % ---- /X is a sparse group: data / indices / indptr --------------------
    idx = find(groupNames == "/X");

    data = pkg.e_guessh5field(filenm, {'/X/'}, {'data'}, true);
    shapeGroupIdx = idx;  % default: read shape from /X
    rawIdx = [];
    rawXIdx = [];
    sparsePath = '/X';    % group whose encoding-type describes the arrays
    probe = data(1:min(5, numel(data)));
    if isequal(probe, round(probe))
        indices = pkg.e_guessh5field(filenm, {'/X/'}, {'indices'}, true);
        indptr = pkg.e_guessh5field(filenm, {'/X/'}, {'indptr'}, true);
    else
        warning('sc_readh5adfile:NormalizedX', ...
            '/X appears transformed/normalized. Attempting to read /raw/X instead.');
        try
            data = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'data'}, true);
            indices = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indices'}, true);
            indptr = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indptr'}, true);
            sparsePath = '/raw/X';
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

    % reconstruct genes-by-cells sparse matrix from the CSR/CSC arrays
    if isMATLABReleaseOlderThan('R2025a')
        data = double(data);   % single-valued sparse needs R2025a or newer
    end
    X = buildSparseFromCS(data, indices, indptr, shape(1), shape(2), ...
        readEncodingType(filenm, sparsePath));
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

% Both are optional metadata and stay empty when absent. The obs table is
% decoded once by PKG.I_READH5ADOBS, and PKG.I_H5ADOBSROLES says which of
% its columns these two outputs come from - the same answer
% PKG.I_ADDH5ADOBSATTRIBS uses to leave them out of the cell attributes, so
% a column cannot land on the object twice under two spellings.
[obsNames, obsValues] = pkg.i_readh5adobs(filenm);
roles = pkg.i_h5adobsroles(obsNames);

if roles.batch > 0
    batchid = obsValues{roles.batch};
end
if roles.celltype > 0
    celltype = obsValues{roles.celltype};
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


function enc = readEncodingType(h5file, groupPath)
% READENCODINGTYPE Read the AnnData "encoding-type" attribute of a group.
%   Returns an empty string when the attribute is absent, as in files written
%   by pre-0.8 AnnData or by h5sparse. Callers must cope with that.

    enc = "";
    try
        enc = lower(strtrim(string(h5readatt(h5file, groupPath, 'encoding-type'))));
    catch
        % Attribute is absent; the caller infers the layout from indptr instead
    end
end


function X = buildSparseFromCS(data, indices, indptr, nCells, nGenes, encoding)
% BUILDSPARSEFROMCS Build a genes-by-cells matrix from CSR or CSC arrays.
%   AnnData stores /X as an (n_obs x n_vars) cells-by-genes matrix, in either
%   CSR or CSC layout. scGEAToolbox works with genes-by-cells, so both layouts
%   are transposed here on the way in.
%
%   The layout is taken from the length of indptr, which is n_obs+1 for CSR and
%   n_vars+1 for CSC. That is decisive except for a square matrix, where the
%   "encoding-type" attribute settles it instead.

    nPtr = numel(indptr) - 1;
    if nCells ~= nGenes
        isCSC = (nPtr == nGenes);
    else
        isCSC = (encoding == "csc_matrix");
    end

    if nPtr ~= nCells && nPtr ~= nGenes
        error('sc_readh5adfile:BadIndptr', ...
            ['indptr has %d entries, matching neither the %d cells nor the ' ...
            '%d genes given by the shape attribute. Check that the file is ' ...
            'complete and stores /X as a CSR or CSC matrix.'], ...
            nPtr, nCells, nGenes);
    end

    % Expand indptr into one major-axis subscript per stored value. AnnData
    % indices are 0-based, so the minor axis shifts up by one.
    major = repelem((1:nPtr)', double(diff(double(indptr(:)))));
    minor = double(indices(:))+1;

    if isCSC
        % Major axis runs over genes, indices point at cells
        X = sparse(major, minor, data(:), nGenes, nCells);
    else
        % Major axis runs over cells, indices point at genes
        X = sparse(minor, major, data(:), nGenes, nCells);
    end
end
