function [X,genelist,keptidxv]=sc_qcfilter(X,genelist,libszcutoff,...
                    mtratio,min_cells_nonzero,gnnumcutoff)
%Basic QC filter
if nargin<6, gnnumcutoff=200; end
if nargin<5 || isempty(min_cells_nonzero)
    min_cells_nonzero=10;         % 0.01
end
% if nargin<5 || isempty(dropout), dropout=1; end           % remove genes no expression
if nargin<4 || isempty(mtratio), mtratio=0.15; end          % 0.10
if nargin<3 || isempty(libszcutoff), libszcutoff=500; end   % 1000

[X,keptidx]=sc_rmmtcells(X,genelist,mtratio);
keptidxv{1}=keptidx;
%if removemtgenes
%    [X,genelist]=sc_rmmtgenes(X,genelist);
%end
oldsz=0;
newsz=1;
c=1;
while ~isequal(oldsz,newsz)
    oldsz=size(X);
    [X,genelist]=sc_filterg(X,genelist);   % remove empty genes
    [X,keptidx]=sc_filterc(X);             % remove empty cells
    keptidxv{end+1}=keptidx;
    [X,genelist]=sc_selectg(X,genelist,min_cells_nonzero);
    [X,genelist]=sc_rmdugenes(X,genelist);
    [X,keptidx]=sc_selectc(X,libszcutoff,gnnumcutoff);
    keptidxv{end+1}=keptidx;
    newsz=size(X);
    c=c+1;
end
end


