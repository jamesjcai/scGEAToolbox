function [X,genelist,keptidxv]=sc_qcfilter(X,genelist,lbsz,mtratio,zeroratio)

if nargin<3, lbsz=1000; end
if nargin<4, mtratio=0.1; end
if nargin<5, zeroratio=0.05; end

[X,keptidx]=sc_rmmtcells(X,genelist,mtratio);
keptidxv{1}=keptidx;
[X,genelist]=sc_rmmtgenes(X,genelist);
oldsz=0;
newsz=1;
c=1;
while ~isequal(oldsz,newsz)
    oldsz=size(X);
    [X,genelist]=sc_filterg(X,genelist,zeroratio);
    [X,keptidx]=sc_filterc(X);
    keptidxv{end+1}=keptidx;
    [X,genelist]=sc_selectg(X,genelist);
    [X,keptidx]=sc_selectc(X,lbsz);
    keptidxv{end+1}=keptidx;
    newsz=size(X);
    c=c+1;
end
