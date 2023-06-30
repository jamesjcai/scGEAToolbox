function [X,g,b,filenm]=sc_readloomfile(filenm)
%Read LOOM file
% http://linnarssonlab.org/loompy/index.html

X=[]; g=[]; b=[];
if nargin<1
[filenm, pathname] = uigetfile( ...
       {'*.loom', 'LOOM Files (*.loom)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a LOOM file');
	if isequal(filenm,0), return; end
	filenm=fullfile(pathname,filenm);
end
if exist(filenm,'file') ~= 2
    error('FileNotFound');
end

fw=gui.gui_waitbar;

%hinfo=h5info(filenm);
% X=h5read(filenm,'/matrix/');
X=pkg.e_guessh5field(filenm,{'/matrix/'},{''},true);

%g=h5read(filenm,'/row_attrs/Gene');
g=pkg.e_guessh5field(filenm,{'/row_attrs/'},{'Gene','gene_ids','gene_name'},false);
if isempty(g), warning('G is not assigned.'); end


b=pkg.e_guessh5field(filenm,{'/col_attrs/'},{'CellID','Cell'},false);
if isempty(b), warning('B is not assigned.'); end

% barcodes=[];
% try
%     barcodes=h5read(filenm,'/col_attrs/Cell');
% catch
%     try
%         barcodes=h5read(filenm,'/col_attrs/CellID');
%     catch ME
%         warning(ME.message);        
%     end
% end

X=pkg.e_uint2sparse(X);

% if ~issparse(X) && ~isa(X,'double')
%     try    
%         X=sparse(double(X));
%     catch ME
%         if (strcmp(ME.identifier,'MATLAB:array:SizeLimitExceeded'))
%             disp('Converting X to sparse.');
%             % tic
%             % S=spalloc(size(X,1),size(X,2),nnz(X));
%             % idx=find(X>0);
%             % S(idx)=X(idx);
%             % toc
%             % X=S;
%             a=floor(size(X)./2);
%             x1=sparse(double(X(1:a(1),1:a(2))));
%             x2=sparse(double(X(a(1)+1:end,1:a(2))));
%             x3=sparse(double(X(1:a(1),a(2)+1:end)));
%             x4=sparse(double(X(a(1)+1:end,a(2)+1:end)));
%             X=[x1 x3; x2 x4];
%         end
%     end
% end

gui.gui_waitbar(fw);
end