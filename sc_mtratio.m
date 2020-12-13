function [r]=sc_mtratio(X,genelist,vebrose)

if nargin<3, vebrose=true; end
r=nan(size(X,2),1);
idx=startsWith(genelist,'mt-','IgnoreCase',true);
if sum(idx)>0
    if vebrose
        fprintf('%d mt-genes found.\n',sum(idx));
    end
    lbsz=sum(X,1);
    lbsz_mt=sum(X(idx,:),1);
    r=transpose(lbsz_mt./lbsz);
else
    if vebrose
        fprintf('No mt-genes found.\n');
    end
end
