function [X,genelist,barcodes,filenm]=sc_readh5adfile(filenm)
%Read HDF5 file
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices

% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3489183
% h5file='GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5';

barcodes=[];
if nargin<1
[filenm, pathname] = uigetfile( ...
       {'*.h5ad', 'H5AD Files (*.h5)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a H5AD file');
	if isequal(filenm,0), X=[]; genelist=[]; return; end
	filenm=fullfile(pathname,filenm);
end
if exist(filenm,'file') ~= 2
    error('FileNotFound');
end

hinfo=h5info(filenm);

idx=find(strcmp(strtrim(string(char(hinfo.Groups.Name))),"/X"));
data=h5read(filenm,[hinfo.Groups(idx).Name,'/data']);
indices=h5read(filenm,[hinfo.Groups(idx).Name,'/indices']);
indptr=h5read(filenm,[hinfo.Groups(idx).Name,'/indptr']);



% idx=find(strcmp(strtrim(string(char(hinfo.Groups.Name))),"/raw"));
% data=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/data']);
% indices=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/indices']);
% indptr=h5read(filenm,[hinfo.Groups(idx).Groups(1).Name,'/indptr']);

idx2=find(strcmp(strtrim(string(char(hinfo.Groups(idx).Attributes.Name))),"shape"));
shape=double(hinfo.Groups(idx).Attributes(idx2).Value);

    idx=find(strcmp(strtrim(string(char(hinfo.Groups.Name))),"/var"));
    g=h5read(filenm,[hinfo.Groups(idx).Name,'/gene_ids']);
    % g=h5read(filenm,[hinfo.Groups(idx).Name,'/gene_name']);

try
    barcodes=h5read(filenm,[hinfo.Groups.Groups(1).Name,'/barcodes']);
catch
        try
barcodes=h5read(filenm,[hinfo.Groups(2).Name,'/index']);
            barcodes=h5read(filenm,[hinfo.Groups(1).Name,'/barcodes']);
        catch
            warning('BARCODES not found.');
        end
end

X=zeros(shape(1),shape(2));
for k=1:length(indptr)-1
    i=indptr(k)+1:indptr(k+1);
    y=indices(i)+1;
    X(k,y)=data(i);
end
X=X.';

genelist=deblank(string(g));
% genelist=strings(length(g),1);
% for k=1:length(g)
%     genelist(k)=string(g(k).data);
% end

end