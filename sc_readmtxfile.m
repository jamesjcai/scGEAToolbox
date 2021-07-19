function [X,genelist,celllist]=sc_readmtxfile(matrixmtxfile,featurestxtfile,barcodestxtfile,coln)
%Read MTX file
if nargin<4, coln=2; end
if nargin<2, featurestxtfile=[]; end

%[X,genelist,celllist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv')
if exist(matrixmtxfile,'file') ~= 2
    error(message('FileNotFound'));        
end
%tic
fprintf('Reading mtx file %s...',matrixmtxfile);
X=pkg.mmread(matrixmtxfile);
try
    X=full(X);
catch

end
fprintf('...done.\n');
if isempty(featurestxtfile)
    genelist=[];
    celllist=[];
    return;
end
fprintf('Reading tsv file %s...',featurestxtfile);
T=readtable(featurestxtfile,'ReadVariableNames',false,...
    'filetype','text','Delimiter',{'\t',',',' ',';','|'});
if coln==1
    genelist=string(T.Var1);
elseif coln==2
    try
        genelist=string(T.Var2);
    catch
        genelist=string(T.Var1);
    end
end
fprintf('...done.\n');
if nargout>2 && ~isempty(barcodestxtfile)
    fprintf('Reading tsv file %s...',barcodestxtfile);
    T=readtable(barcodestxtfile,'FileType','text','Delimiter','\n','ReadVariableNames',false);
    celllist=string(T.Var1);
    fprintf('...done.\n');
end
%assert(isequal(size(X,1),length(genelist)))
%toc
end