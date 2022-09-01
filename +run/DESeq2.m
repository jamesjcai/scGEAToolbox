function [T,Tup,Tdn]=DESeq2(X,Y,gene)

if nargin<3, gene=(1:size(X,1))'; end
pw1=fileparts(mfilename('fullpath'));
pth=fullfile(pw1,'thirdparty','DESeq2');
if ~(ismcc || isdeployed)
    addpath(pth);
end

avg_1 = mean(X,2);
avg_2 = mean(Y,2);
pct_1 = sum(X>0,2)./size(X,2);
pct_2 = sum(Y>0,2)./size(Y,2);


[avg_log2FC,p_val_adj,~,p_val] = DESeq2ori(X',Y');
abs_log2FC=abs(avg_log2FC);
% T = table(gene, pvalue, abs(log2FC), log2FC, FDR, 'VariableNames', ...
%     ["gene", "p_val", "abs_log2FC", "log2FC", "p_val_adj"]);
% T=sortrows(T,'abs_log2FC','descend');
% T=sortrows(T,'p_val_adj','ascend');

    T = table(gene, p_val, avg_log2FC, abs_log2FC, avg_1, avg_2, ...
           pct_1, pct_2, p_val_adj);

    if nargout>1
        [Tup,Tdn]=pkg.e_processDETable(T);
    end
    
end