function [X,genelist,idx]=sc_rmmtgenes(X,genelist,mtgenenamepat,vebrose)
if nargin<3
   mtgenenamepat="mt-";
end
if nargin<4, vebrose=true; end
idx=startsWith(genelist,mtgenenamepat,'IgnoreCase',true);
if vebrose
    genelist(idx);
end
genelist=genelist(~idx);
X=X(~idx,:);
if sum(idx)>0
   if vebrose
   fprintf('%d mt-genes found and removed.\n',sum(idx));
   end
else
   if vebrose
   fprintf('No mt-genes found or removed.\n',sum(idx));
   end
end
if nargout>3
    idx=find(idx);
end

% find(contains(genelist,"ENSG00000198804"))  % MT-CO1