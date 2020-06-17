function [X,keptidx]=sc_rmmtcells(X,genelist,mtratio,txtpat,vebrose)

if nargin<3, mtratio=0.1; end
if nargin<4, txtpat="mt-"; end
if nargin<5, vebrose=true; end

idx=startsWith(genelist,txtpat,'IgnoreCase',true);
if sum(idx)>0 && vebrose
   fprintf('%d mt-genes found.\n',sum(idx));
else
   fprintf('No mt-genes found.\n');
end
lbsz=sum(X); 
lbsz_mt=sum(X(idx,:));
f_mtreads=lbsz_mt./lbsz;
keptidx=f_mtreads<mtratio;
if sum(~keptidx)>0
    X=X(:,keptidx);
    if vebrose && ~issparse(X)
        fprintf('%d cells with >=%.2f mt-read ratio removed.\n',...
            sum(~keptidx),mtratio);
    end
end
