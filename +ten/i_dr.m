function [T] = i_dr(aln0, aln1, genelist, dosort, reference)
% DR - differential regulatory gene identification
%
% T = ten.i_dr(aln0, aln1, genelist, dosort, reference)
%
% REFERENCE (default false) follows dRegulation() of the R packages. The
% p-values are the same either way (FC = d^2/mean(d^2) in both); reference
% mode adds what R has on top of them:
%   distance - the Euclidean distance d itself (drdist is d^2/norm(d^2))
%   Z        - d raised to the Box-Cox power MASS::boxcox selects on
%              seq(-2, 2, length.out = 1000), then standardized
% sets pValues to 1 for genes whose distance is at the level of
% floating-point noise (d <= sqrt(eps)*max|alignment|), as R does, and sorts
% by pValues ascending, keeping ties in gene order, instead of by drdist.
if nargin < 5, reference = false; end
if nargin < 4, dosort = true; end
if nargin < 3, genelist = string(num2cell(1:size(aln0, 1)))'; end
drdist = vecnorm(aln0-aln1, 2, 2).^2;
drdist = drdist ./ norm(drdist);
FC = drdist ./ mean(drdist);
pValues = chi2cdf(FC, 1, 'upper');
if reference
    distance = vecnorm(aln0-aln1, 2, 2);
    Z = i_boxcoxz(distance);
    isNoise = distance <= sqrt(eps)*max(abs([aln0; aln1]), [], "all");
    pValues(isNoise) = 1;
    if all(isNoise)
        warning('ten:i_dr:allNoise', ...
            ['No gene differs between the two conditions beyond numerical ' ...
             'noise; all p-values were set to 1.']);
    end
end

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
if reference
    T = addvars(T, distance, Z, 'After', 'genelist');
    if dosort
        T = sortrows(T, 'pValues', 'ascend');
    end
elseif dosort
    T = sortrows(T, 'drdist', 'descend');
end
end

function Z = i_boxcoxz(d)
% scale(nD) of dRegulation, nD being d at the Box-Cox power of largest
% profile log-likelihood (MASS::boxcox(d ~ 1)). If a distance is not
% positive boxcox fails, and R standardizes d itself.
nD = d;
if all(d > 0)
    lambdas = linspace(-2, 2, 1000);
    lambdas = lambdas(lambdas ~= 0);
    y = d/exp(mean(log(d)));    % geometric-mean scaling, as MASS
    logy = log(y);
    n = numel(y);
    logLik = zeros(size(lambdas));
    for k = 1:numel(lambdas)
        la = lambdas(k);
        if abs(la) > 1/50
            yt = (y.^la - 1)/la;
        else
            yt = logy.*(1 + (la*logy)/2.*(1 + (la*logy)/3.*(1 + (la*logy)/4)));
        end
        logLik(k) = -n/2*log(sum((yt - mean(yt)).^2));
    end
    [~, best] = max(logLik);
    bc = lambdas(best);
    if bc < 0
        nD = 1./(d.^bc);
    else
        nD = d.^bc;
    end
end
Z = (nD - mean(nD))/std(nD);
end
