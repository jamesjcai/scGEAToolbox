function [X]=sc_rmmtcells(X,genelist,mtratio,txtpat,vebrose)

if nargin<3, mtratio=0.1; end
if nargin<4, txtpat="mt-"; end
if nargin<5, vebrose=true; end

idx=startsWith(genelist,txtpat,'IgnoreCase',true);
if sum(idx)>0 && vebrose
   fprintf('%d mt-genes found.\n',sum(idx));
else
   fprintf('No mt-genes found.\n',sum(idx));  
end
lbsz=sum(X); 
lbsz_mt=sum(X(idx,:));
f_mtreads=lbsz_mt./lbsz;
if sum(f_mtreads>=mtratio)>0
    X=X(:,f_mtreads<mtratio);
    if vebrose
        fprintf('%d cells with >=%.2f mt-read ratio removed.\n',...
            sum(f_mtreads<mtratio),mtratio);
    end
end
