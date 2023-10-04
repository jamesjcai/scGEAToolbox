function [X, g, b, filenm] = sc_readh5adfile(filenm)
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

if nargin < 1 || isempty(filenm)
    [filenm, pathname] = uigetfile( ...
        {'*.h5ad', 'H5AD Files (*.h5ad)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a H5AD file');
    if isequal(filenm, 0), return; end
    filenm = fullfile(pathname, filenm);
end
if exist(filenm, 'file') ~= 2, error('File Not Found.'); end

fw = gui.gui_waitbar_adv;

hinfo = h5info(filenm);

idx = find(strcmp(strtrim(string(char(hinfo.Groups.Name))), "/X"));
%data=h5read(filenm,[hinfo.Groups(idx).Name,'/data']);
%indices=h5read(filenm,[hinfo.Groups(idx).Name,'/indices']);
%indptr=h5read(filenm,[hinfo.Groups(idx).Name,'/indptr']);

data = pkg.e_guessh5field(filenm, {'/X/'}, {'data'}, true);
if isequal(data(1:5), round(data(1:5)))
    indices = pkg.e_guessh5field(filenm, {'/X/'}, {'indices'}, true);
    indptr = pkg.e_guessh5field(filenm, {'/X/'}, {'indptr'}, true);
else
    disp('Reading /raw/X');
    data = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'data'}, true);
    indices = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indices'}, true);
    indptr = pkg.e_guessh5field(filenm, {'/raw/X/'}, {'indptr'}, true);
end


% idx=find(strcmp(strtrim(string(char(hinfo.Groups.Name))),"/raw"));
% data=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/data']);
% indices=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/indices']);
% indptr=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/indptr']);

idx2 = find(strcmp(strtrim(string(char(hinfo.Groups(idx).Attributes.Name))), "shape"));
shape = double(hinfo.Groups(idx).Attributes(idx2).Value);

g = pkg.e_guessh5field(filenm, {'/var/'}, {'_index', 'gene_ids', 'gene_name'}, false);

if isempty(g) || length(unique(strlength(g))) == 1 % suggesting ENSEMBLE ID
    disp('Reading /var/feature_name/categories');
    gx = pkg.e_guessh5field(filenm, {'/raw/var/feature_name/', ...
        '/var/feature_name/'}, {'categories'}, false);
    if ~isempty(gx)
        g = gx;
    end
end
if isempty(g), warning('G is not assigned.'); end


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

b = pkg.e_guessh5field(filenm, {'/obs/'}, {'_index', 'barcodes','cell_id'});
if isempty(b), warning('B is not assigned.'); end


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

X = spalloc(shape(2), shape(1), length(data));
for k = 1:length(indptr) - 1
    i = indptr(k) + 1:indptr(k+1);
    y = indices(i) + 1;
    X(y, k) = data(i);
end

g = deblank(string(g));

% genelist=strings(length(g),1);
% for k=1:length(g)
%     genelist(k)=string(g(k).data);
% end
gui.gui_waitbar_adv(fw);

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