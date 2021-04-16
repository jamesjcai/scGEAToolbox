function [X,genelist,keptidxv]=sc_qcfilter(X,genelist,libsize,mtratio,dropout,min_cells_nonzero,removemtgenes)

if nargin<7, removemtgenes=false; end
if nargin<6 || isempty(min_cells_nonzero), min_cells_nonzero=0.01; end
if nargin<5 || isempty(dropout), dropout=0.01; end
if nargin<4 || isempty(mtratio), mtratio=0.1; end
if nargin<3 || isempty(libsize), libsize=1000; end

[X,keptidx]=sc_rmmtcells(X,genelist,mtratio);
keptidxv{1}=keptidx;
if removemtgenes
    [X,genelist]=sc_rmmtgenes(X,genelist);
end
oldsz=0;
newsz=1;
c=1;
while ~isequal(oldsz,newsz)
    oldsz=size(X);
    [X,genelist]=sc_filterg(X,genelist,dropout);
    [X,keptidx]=sc_filterc(X);
    keptidxv{end+1}=keptidx;
    [X,genelist]=sc_selectg(X,genelist,1,min_cells_nonzero);
    [X,genelist]=sc_rmdugenes(X,genelist);
    [X,keptidx]=sc_selectc(X,libsize);
    keptidxv{end+1}=keptidx;
    newsz=size(X);
    c=c+1;
end
end


