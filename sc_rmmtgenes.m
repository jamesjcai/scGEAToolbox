function [X,genelist,idx]=sc_rmmtgenes(X,genelist,txtpat,vebrose)
if nargin<3
   txtpat="mt-";
end
if nargin<4, vebrose=true; end
idx=startsWith(genelist,txtpat,'IgnoreCase',true);
% genelist(idx);
genelist=genelist(~idx);
X=X(~idx,:);
if sum(idx)>0 && vebrose
   fprintf('%d mt-genes removed.\n',sum(idx));
else
   fprintf('No mt-genes found.\n',sum(idx));
end
if nargout>3
    idx=find(idx);
end

% find(contains(genelist,"ENSG00000198804"))  % MT-CO1