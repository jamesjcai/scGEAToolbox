function [X,genelist]=sc_readh5file(h5filename)

% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices

% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3489183
% h5file='GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5';
if nargin<1
[h5filename, pathname] = uigetfile( ...
       {'*.h5;*.hdf5', 'HDF5 Files (*.h5)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a HDF5 file');
	if isequal(h5filename,0), X=[]; genelist=[]; return; end
	h5filename=fullfile(pathname,h5filename);
end
if exist(h5filename,'file') ~= 2
    error(message('FileNotFound'));        
end

% h5disp(h5file);
hinfo=h5info(h5filename);

% if strcmp(a.Groups(1).Datasets(2).Name,'data')
data=h5read(h5filename,[hinfo.Groups(1).Name,'/data']);
indices=h5read(h5filename,[hinfo.Groups(1).Name,'/indices']);
indptr=h5read(h5filename,[hinfo.Groups(1).Name,'/indptr']);
g=h5read(h5filename,[hinfo.Groups(1).Name,'/gene_names']);
shape=h5read(h5filename,[hinfo.Groups(1).Name,'/shape']);

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