function [X,genelist]=sc_qcfilter(X,genelist,libsizecutoff)

if nargin<3
    libsizecutoff=20000;
end
[X,genelist]=sc_rmmtgenes(X,genelist);
oldsz=0;
newsz=1;
c=0;
while ~isequal(oldsz,newsz)
    oldsz=size(X);
    [X,genelist]=sc_filterg(X,genelist);
    [X]=sc_filterc(X);
    [X,genelist]=sc_selectg(X,genelist);
    [X]=sc_selectc(X,libsizecutoff);
    newsz=size(X);
    c=c+1;
end
c

