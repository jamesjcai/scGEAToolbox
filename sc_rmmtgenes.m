function [X,genelist]=sc_rmmtgenes(X,genelist)

idx=startsWith(genelist,'mt-','IgnoreCase',true);
genelist(idx)
genelist=genelist(~idx);
X=X(~idx,:);
