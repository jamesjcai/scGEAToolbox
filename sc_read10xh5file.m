function [X, g, b, c, filenm] = sc_read10xh5file(filenm)
%Read 10x Genomics H5 file
% https://www.mathworks.com/help/matlab/hdf5-files.html
% http://scipy-lectures.org/advanced/scipy_sparse/csc_matrix.html
% https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/h5_matrices

% https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3489183
% h5file='GSM3489183_IPF_01_filtered_gene_bc_matrices_h5.h5';

X = []; % expression matrix
g = []; % gene names
b = []; % barcode of cells
c = []; % batch id
if nargin < 1 || isempty(filenm)
    [filenm, pathname] = uigetfile( ...
        {'*.h5;*.hdf5', 'HDF5 Files (*.h5)'; ...
        '*.*', 'All Files (*.*)'}, ...
        'Pick a 10x Genomics H5 file');
    if isequal(filenm, 0), return; end
    filenm = fullfile(pathname, filenm);
end
if exist(filenm, 'file') ~= 2
    error('FileNotFound');
end

grouptag = "/matrix/";
%try
    h = h5info(filenm);
    % grouptag = strcat(h.Groups(1).Name, "/");
    assert(any(contains(string({h.Groups.Name}), "/matrix")))
%catch
%    warning('Failed to read HDF5 file info. Using default group tag.');
%end

data = pkg.e_guessh5field(filenm, {grouptag}, {'data'}, true);
indices = pkg.e_guessh5field(filenm, {grouptag}, {'indices'}, true);
indptr = pkg.e_guessh5field(filenm, {grouptag}, {'indptr'}, true);
shape = pkg.e_guessh5field(filenm, {grouptag}, {'shape'}, true);


g = pkg.e_guessh5field(filenm, {grouptag, '/matrix/features/'}, {'gene_names', 'name'}, false);
if isempty(g)
    warning('Gene names or feature names are not assigned.');
end

% try
% g=h5read(filenm,[h.Groups.Groups(1).Name,'/gene_names']);
% catch
%     try
%         g=h5read(filenm,[h.Groups.Groups(1).Name,'/name']);
%     catch
%         try
%             g=h5read(filenm,[h.Groups(1).Name,'/gene_names']);
%         catch
%             error('GENE_NAMES not found.');
%         end
%     end
% end

b = pkg.e_guessh5field(filenm, {grouptag, '/matrix/features/'}, {'barcodes'}, false);
if isempty(b), warning('B is not assigned.'); end


%
%     try
%     barcodes=h5read(filenm,[hinfo.Groups.Groups(1).Name,'/barcodes']);
% catch
%         try
%             barcodes=h5read(filenm,[hinfo.Groups(1).Name,'/barcodes']);
%         catch
%             warning('BARCODES not found.');
%         end
% end

% try
%     X=zeros(shape(1),shape(2));
% catch
if ~isMATLABReleaseOlderThan('R2025a')
    X = spalloc(shape(1), shape(2), length(data), 'single');
else
    X = spalloc(shape(1), shape(2), length(data));
end

%end

%c=0; olda=-1;
for k = 1:length(indptr) - 1
    % if mod(c,round(length(indptr)/100))==0
    %     a=round(100*(c/length(indptr)));
    %     if a~=olda
    %         %fprintf('......%d%%\n',a);
    %         gui.gui_waitbar_adv(fw,a/100);
    %         olda=a;
    %     end
    % end
    ix = indptr(k) + 1:indptr(k+1);
    X((indices(ix) + 1), k) = data(ix);
    %    c=c+1;
end
%fprintf('......100%%\n');

g = deblank(string(g));

% genelist=strings(length(g),1);
% for k=1:length(g)
%     genelist(k)=string(g(k).data);
% end
%gui.gui_waitbar(fw);
%gui.gui_waitbar_adv(fw);

if all(contains(b,'-'))
    try
        % c = extractAfter(b, 25);
        c = extractAfter(b, cell2mat(strfind(b,"-")));
    catch
    end
end

end



%{

function countMatrix = getMatrixFromH5(filename)
    info = h5info(filename, '/matrix');
    
    barcodes = h5read(filename, '/matrix/barcodes');
    data = h5read(filename, '/matrix/data');
    indices = h5read(filename, '/matrix/indices');
    indptr = h5read(filename, '/matrix/indptr');
    shape = h5read(filename, '/matrix/shape');
    
    matrix = sparse(indices+1, indptr+1, data, shape(2), shape(1));
    
    featureRef = struct();
    featureGroup = info.Groups(strcmp({info.Groups.Name}, '/matrix/features'));
    featureRef.id = h5read(filename, '/matrix/features/id');
    featureRef.name = h5read(filename, '/matrix/features/name');
    featureRef.featureType = h5read(filename, '/matrix/features/feature_type');
    
    tagKeys = h5read(filename, '/matrix/features/_all_tag_keys');
    for i = 1:length(tagKeys)
        key = char(tagKeys(i));
        featureRef.(key) = h5read(filename, ['/matrix/features/' key]);
    end
    
    countMatrix = struct('featureRef', featureRef, 'barcodes', barcodes, 'matrix', matrix);
end

filteredH5 = '/opt/sample345/outs/filtered_feature_bc_matrix.h5';
filteredMatrix = getMatrixFromH5(filteredH5);

% https://www.10xgenomics.com/support/software/space-ranger/advanced/hdf5-feature-barcode-matrix-format

%}


