function [s] = e_fgsearun(T, rmribo)
% Run fast GSEA (fGSEA) analysis in R
if nargin < 2, rmribo = true; end
[s] = run.r_fgsea(T.genelist, rmribo);
end
