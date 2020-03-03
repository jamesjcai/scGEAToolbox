function [X,genelist]=sc_qcfilter(X,genelist,lbsz,mtratio,zeroratio)

if nargin<3, lbsz=20000; end
if nargin<4, mtratio=0.1; end
if nargin<5, zeroratio=0.05; end

[X]=sc_rmmtcells(X,genelist,mtratio);
[X,genelist]=sc_rmmtgenes(X,genelist);
oldsz=0;
newsz=1;
c=0;
while ~isequal(oldsz,newsz)
    oldsz=size(X);
    [X,genelist]=sc_filterg(X,genelist,zeroratio);
    [X]=sc_filterc(X);
    [X,genelist]=sc_selectg(X,genelist);
    [X]=sc_selectc(X,lbsz);
    newsz=size(X);
    c=c+1;
end
