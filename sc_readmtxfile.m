function [X,genelist,celllist]=sc_readmtxfile(matrixmtxfile,featurestxtfile,barcodestxtfile,coln)
if nargin<4, coln=1; end
if nargin<2, featurestxtfile=[]; end

%[X,genelist,celllist]=sc_readmtxfile('matrix.mtx','features.tsv','barcodes.tsv')
if exist(matrixmtxfile,'file') ~= 2
    error(message('FileNotFound'));        
end
X=pkg.mmread(matrixmtxfile);
try
    X=full(X);
catch

end
if isempty(featurestxtfile)
    genelist=[];
    celllist=[];
    return;
end
T=readtable(featurestxtfile,'ReadVariableNames',false,...
    'filetype','text','Delimiter',{'comma','space','tab','semi','bar'});
if coln==1
    genelist=string(T.Var1);
elseif coln==2
    try
        genelist=string(T.Var2);
    catch
        genelist=string(T.Var1);
    end
end
if nargout>2 && ~isempty(barcodestxtfile)
    T=readtable(barcodestxtfile,'ReadVariableNames',false,'filetype','text');
    celllist=string(T.Var1);
end

%assert(isequal(size(X,1),length(genelist)))

end