function [X,genelist]=sc_mmread(mtxfile,genefile)
if nargin<1
[mtxfile, pathname] = uigetfile( ...
       {'*.csv;*.tsv;*.tab;*.txt', 'MatrixMarket Files (*.mtx)';
        '*.*',  'All Files (*.*)'}, ...
        'Pick a MatrixMarket file');
	if isequal(mtxfile,0), X=[]; genelist=[]; return; end
	mtxfile=fullfile(pathname,mtxfile);
end
if exist(mtxfile,'file') ~= 2
    error(message('FileNotFound'));        
end
if nargin>1
    if exist(genefile,'file') ~= 2
        error(message('FileNotFound'));
    end
end
[X,nrows]=mmread(mtxfile);
genelist=[];

if nargin>1
    T=readtable(genefile,'HeaderLines',0,'filetype','text',...
                'ReadVariableNames',false,'Delimiter','\t,',...
                'Format','%s');
    genelist=string(table2array(T(:,1)));
    assert(nrows==length(genelist));
end
