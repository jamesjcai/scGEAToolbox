function [X,genelist,idx]=sc_rmmtgenes(X,genelist,txtpat)
if nargin<3
   txtpat="mt-";
end
idx=startsWith(genelist,txtpat,'IgnoreCase',true);
% genelist(idx);
genelist=genelist(~idx);
X=X(~idx,:);
if nargout>3
    idx=find(idx);
end
