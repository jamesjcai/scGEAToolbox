function [X,genelist,barcodes,filenm]=sc_readloomfile(filenm)
%Read LOOM file
% http://linnarssonlab.org/loompy/index.html

barcodes=[];
if nargin<1
[filenm, pathname] = uigetfile( ...
       {'*.loom', 'LOOM Files (*.loom)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a LOOM file');
	if isequal(filenm,0), X=[]; genelist=[]; return; end
	filenm=fullfile(pathname,filenm);
end
if exist(filenm,'file') ~= 2
    error('FileNotFound');
end

%hinfo=h5info(filenm);

X=h5read(filenm,'/matrix');
genelist=h5read(filenm,'/row_attrs/Gene');

barcodes=[];
try
    barcodes=h5read(filenm,'/col_attrs/Cell');
catch
    try
        barcodes=h5read(filenm,'/col_attrs/CellID');
    catch ME
        warning(ME.message);        
    end
end


try
    
    %X=sparse(double(X));
catch
    % tic
    % S=spalloc(size(X,1),size(X,2),nnz(X));
    % idx=find(X>0);
    % S(idx)=X(idx);
    % toc
    % X=S;
end
end