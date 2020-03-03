function [X]=sc_rmmtcells(X,genelist,mtratio,txtpat)

if nargin<3, mtratio=0.1; end
if nargin<4, txtpat="mt-"; end
idx=startsWith(genelist,txtpat,'IgnoreCase',true);
lbsz=sum(X); 
lbsz_mt=sum(X(idx,:));
f_mtreads=lbsz_mt./lbsz;
X=X(:,f_mtreads<mtratio);
