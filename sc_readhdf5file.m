function [X,genelist,barcodes,filenm]=sc_readhdf5file(filenm)
%Read HDF5 file
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices

% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3489183
% h5file='GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5';

barcodes=[];
if nargin<1
[filenm, pathname] = uigetfile( ...
       {'*.h5;*.hdf5', 'HDF5 Files (*.h5)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a HDF5 file');
	if isequal(filenm,0), X=[]; genelist=[]; return; end
	filenm=fullfile(pathname,filenm);
end
if exist(filenm,'file') ~= 2
    error('FileNotFound');
end

hinfo=h5info(filenm);
% h5disp(filenm,'/matrix','min');

% if strcmp(a.Groups(1).Datasets(2).Name,'data')
data=h5read(filenm,[hinfo.Groups(1).Name,'/data']);
indices=h5read(filenm,[hinfo.Groups(1).Name,'/indices']);
indptr=h5read(filenm,[hinfo.Groups(1).Name,'/indptr']);
shape=h5read(filenm,[hinfo.Groups(1).Name,'/shape']);

try
g=h5read(filenm,[hinfo.Groups.Groups(1).Name,'/gene_names']);
catch
    try
        g=h5read(filenm,[hinfo.Groups.Groups(1).Name,'/name']);
    catch
        try
            g=h5read(filenm,[hinfo.Groups(1).Name,'/gene_names']);
        catch
            error('GENE_NAMES not found.');
        end
    end
end

try
    barcodes=h5read(filenm,[hinfo.Groups.Groups(1).Name,'/barcodes']);
catch
        try
            barcodes=h5read(filenm,[hinfo.Groups(1).Name,'/barcodes']);
        catch
            warning('BARCODES not found.');
        end
end

X=zeros(shape(1),shape(2));
for k=1:length(indptr)-1
    i=indptr(k)+1:indptr(k+1);
    y=indices(i)+1;
    X(y,k)=data(i);
end


genelist=deblank(string(g));
% genelist=strings(length(g),1);
% for k=1:length(g)
%     genelist(k)=string(g(k).data);
% end

end