function [X, g, b, batchid, celltype, filenm] = sc_readh5adfile(filenm)
%Read H5AD file
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


%names = string(pkg.i_extractfield(hinfo.Groups, 'Name'));
%idx = find(names == "/X");

% idx = find(string({hinfo.Groups.Name})=="/X")

%data=h5read(filenm,[hinfo.Groups(idx).Name,'/data']);
%indices=h5read(filenm,[hinfo.Groups(idx).Name,'/indices']);
%indptr=h5read(filenm,[hinfo.Groups(idx).Name,'/indptr']);

data = pkg.e_guessh5field(filenm, {'/X/'}, {'data'}, true);
if isequal(data(1:5), round(data(1:5)))
    indices = pkg.e_guessh5field(filenm, {'/X/'}, {'indices'}, true);
    indptr = pkg.e_guessh5field(filenm, {'/X/'}, {'indptr'}, true);
else
    disp('/X has been transformed/normalized. Try to read /raw/X...');
    try
        data = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'data'}, true);
        indices = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indices'}, true);
        indptr = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indptr'}, true);
        disp('/raw/X is read.');
    catch ME
        disp(ME.message)
        indices = pkg.e_guessh5field(filenm, {'/X/'}, {'indices'}, true);
        indptr = pkg.e_guessh5field(filenm, {'/X/'}, {'indptr'}, true);
        disp('Reading transformed/normalized /X instead.');
    end
end


% idx=find(strcmp(strtrim(string(char(hinfo.Groups.Name))),"/raw"));
% data=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/data']);
% indices=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/indices']);
% indptr=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/indptr']);

idx2 = find(strcmp(strtrim(string(char(hinfo.Groups(idx).Attributes.Name))), "shape"));
if isempty(idx2)
idx2 = find(strcmp(strtrim(string(char(hinfo.Groups(idx).Attributes.Name))), "h5sparse_shape"));
end
shape = double(hinfo.Groups(idx).Attributes(idx2).Value);

g = pkg.e_guessh5field(filenm, {'/var/'}, {'_index', 'gene_ids', ...
    'gene_name','symbol'}, false);

if isempty(g) || isscalar(unique(strlength(g))) % suggesting ENSEMBLE ID
    disp('Reading /var/feature_name/categories');
    gx = pkg.e_guessh5field(filenm, {'/raw/var/feature_name/', ...
        '/var/feature_name/'}, {'categories'}, false);
    if ~isempty(gx)
        g = gx;
    end
end
if isempty(g), warning('Genename is not assigned.'); end


% idx=find(strcmp(strtrim(string(char(hinfo.Groups.Name))),"/var"));
% try
%     g=h5read(filenm,[hinfo.Groups(idx).Name,'/_index']);
% catch
%     try
%         g=h5read(filenm,[hinfo.Groups(idx).Name,'/gene_ids']);
%     catch
%         try
%             g=h5read(filenm,[hinfo.Groups(idx).Name,'/gene_name']);
%         catch
%             warning('GENELIST not found.');
%         end
%     end
% end

b = pkg.e_guessh5field(filenm, {'/obs/'}, {'_index', 'barcodes','cell_id','CellID'});
if isempty(b), warning('Barcode is not assigned.'); end

try
    % c = pkg.e_guessh5field(filenm, {'/obs/BatchID/'}, {'codes'});
    % cL = pkg.e_guessh5field(filenm, {'/obs/BatchID/'}, {'categories'});
    % batchid = cL(c+1);
    batchid = readObsColumn(filenm, 'BatchID');
catch ME
    disp(ME.message)
end

try
    % c = pkg.e_guessh5field(filenm, {'/obs/CellType/'}, {'codes'});
    % cL = pkg.e_guessh5field(filenm, {'/obs/CellType/'}, {'categories'});
    % celltype = cL(c+1);
    celltype = readObsColumn(filenm, 'CellType');
catch ME
    disp(ME.message)
end

% try
%     barcodes=h5read(filenm,[hinfo.Groups.Groups(1).Name,'/barcodes']);
% catch
%         try
%             barcodes=h5read(filenm,[hinfo.Groups(2).Name,'/index']);
%             %barcodes=h5read(filenm,[hinfo.Groups(1).Name,'/barcodes']);
%         catch
%             try
%                 barcodes=h5read(filenm,'/obs/_index');
%             catch
%                 warning('BARCODES not found.');
%             end
%         end
% end


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

% genelist=strings(length(g),1);
% for k=1:length(g)
%     genelist(k)=string(g(k).data);
% end
% gui.gui_waitbar_adv(fw);

end


% a=h5read(filenm,'/obs/seurat_clusters')
% a=h5read(filenm,'/var/_index')
% a=h5read(filenm,'/obs/annotation')

%{
function i_tryread
try
    g=h5read(filenm,[hinfo.Groups(idx).Name,'/_index']);
catch

end
end
%}


function values = readObsColumn(h5file, colname)
%READOBSCOLUMN Extract a column from AnnData .h5ad obs table.
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
        if ischar(cats)
            cats = string(cats);
        elseif iscellstr(cats)
            cats = string(cats);
        elseif isstring(cats)
            % ___
        end
        
        % AnnData categorical codes are 0-based, MATLAB is 1-based
        values = categorical(codes + 1, 1:numel(cats), cats);
        return
    catch
        %error('Column "%s" not found in %s\n(HDF5 error: %s)', ...
        %      colname, h5file, ME.message);
    end
end
