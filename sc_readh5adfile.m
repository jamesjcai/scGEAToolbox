function [X, g, b, batchid, celltype, filenm] = sc_readh5adfile(filenm)
% Read H5AD file
% https://anndata.readthedocs.io/en/latest/fileformat-prose.html
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices

% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3489183
% h5file='GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5';

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

% fw = gui.gui_waitbar_adv;

hinfo = h5info(filenm);

idx = find(strcmp(strtrim(string(char(hinfo.Groups.Name))), "/X"));


% names = string(pkg.i_extractfield(hinfo.Groups, 'Name'));
% idx = find(names == "/X");

% idx = find(string({hinfo.Groups.Name})=="/X")

% data=h5read(filenm,[hinfo.Groups(idx).Name,'/data']);
% indices=h5read(filenm,[hinfo.Groups(idx).Name,'/indices']);
% indptr=h5read(filenm,[hinfo.Groups(idx).Name,'/indptr']);

data = pkg.e_guessh5field(filenm, {'/X/'}, {'data'}, true);
shapeGroupIdx = idx;  % default: read shape from /X
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

g = pkg.e_guessh5field(filenm, {'/var/'}, {'_index', 'gene_ids', ...
    'gene_name','symbol'}, false);

if isempty(g) || isscalar(unique(strlength(g))) % suggesting ENSEMBLE ID
    % Gene IDs look uniform-length (e.g. ENSEMBL); try feature_name for symbols
    gx = pkg.e_guessh5field(filenm, {'/raw/var/feature_name/', ...
        '/var/feature_name/'}, {'categories'}, false);
    if ~isempty(gx)
        g = gx;
    end
end
if isempty(g), warning('sc_readh5adfile:NoGenenames', 'Genename is not assigned.'); end

b = pkg.e_guessh5field(filenm, {'/obs/'}, {'_index', 'barcodes','cell_id','CellID'});
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

g = deblank(string(g));

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
