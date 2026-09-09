function [T] = i_dr(aln0, aln1, genelist, dosort)
% DR - differential regulatory gene identification
if nargin < 4, dosort = true; end
if nargin < 3, genelist = string(num2cell(1:size(aln0, 1)))'; end
drdist = vecnorm(aln0-aln1, 2, 2).^2;
drdist = drdist ./ norm(drdist);
FC = drdist ./ mean(drdist);
pValues = chi2cdf(FC, 1, 'upper');

% PKG.E_FDR replaces the two-branch block that used to sit here. Its MAFDR
% branch and its PKG.E_FDR_BH branch are the same algorithm on clean input
% -- they agree to 2e-16 -- but they disagree whenever the p-values carry
% NaN, which is what RANKSUM returns for a gene with no counts in either
% group and what a per-cell-type run produces in bulk. One drops those from
% the family, the other counts them, so the same command gave a different
% answer depending on whether the Bioinformatics Toolbox was installed.
pAdjusted = pkg.e_fdr(pValues);
% if size(genelist,1)==1, genelist=genelist'; end
genelist = genelist(:);
sortid = (1:length(genelist))';
if size(genelist, 2) > 1, genelist = genelist'; end
T = table(sortid, genelist, drdist, FC, pValues, pAdjusted);
if dosort
    T = sortrows(T, 'drdist', 'descend');
end
end
