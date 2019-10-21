function [X,genelist,celllist]=sc_readmtxfile(matrixmtxfile,featurestxtfile,barcodestxtfile)
%[X,genelist,celllist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv')
if exist(matrixmtxfile,'file') ~= 2
    error(message('FileNotFound'));        
end
X=mmread(matrixmtxfile);
T=readtable(featurestxtfile,'ReadVariableNames',false,'filetype','text');
genelist=string(T.Var1);
if nargout>2
    T=readtable(barcodestxtfile,'ReadVariableNames',false,'filetype','text');
    celllist=string(T.Var1);
end
X=full(X);
assert(isequal(size(X,1),length(genelist)))
assert(isequal(size(X,1),length(genelist)))

