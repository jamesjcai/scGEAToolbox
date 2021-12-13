function [T,Tup,Tdn] = sc_deg(X, Y, genelist, methodid)
%DEG analysis using Mann–Whitney U test
    % https://satijalab.org/seurat/v3.1/de_vignette.html
    % p_val : p_val (unadjusted)
    % avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
    % abs_logFC
    % pct.1 : The percentage of cells where the feature is detected in the first group
    % pct.2 : The percentage of cells where the feature is detected in the second group
    % p_val_adj : Adjusted p-value, based on bonferroni correction using all features in the dataset.
    %
    % SEE ALSO: [T]=run.MAST(X,Y,genelist);

    if nargin < 2, error("USAGE: sc_deg(X,Y)\n"); end
    if nargin < 3, genelist = string(1:size(X,1))'; end
    if nargin < 4, methodid = 1; end

    ng = size(X, 1);
    assert(isequal(ng, size(Y, 1)));

    p_val = ones(ng, 1);
    avg_log2FC = ones(ng, 1);
    pct_1 = ones(ng, 1);
    pct_2 = ones(ng, 1);
    
    nx=size(X,2);
    ny=size(Y,2);
    
    X=log(1+sc_norm(X));
    Y=log(1+sc_norm(Y));
    
    for k = 1:ng
        x = X(k, :);
        y = Y(k, :);
        switch methodid
            case 1
                % “wilcox” : Wilcoxon rank sum test (default)
                % also called Mann–Whitney U test
                p_val(k) = ranksum(x, y);
            case 2
                [~, p] = ttest2(x, y);
                p_val(k) = p;
            otherwise
                error('Unknown option')
        end
        avg_log2FC(k) = log2(mean(x) ./ mean(y));
        pct_1(k) = sum(x > 0) ./ nx;
        pct_2(k) = sum(y > 0) ./ ny;
    end

    if exist('mafdr.m', 'file')
        p_val_adj = mafdr(p_val, 'BHFDR', true);
    else
        [~, ~, ~, p_val_adj] = pkg.fdr_bh(p_val);
    end

    %   sortid=(1:length(genelist))';
    if size(genelist, 2) > 1
        gene = genelist';
    else
        gene = genelist;
    end
    abs_log2FC = abs(avg_log2FC);
    T = table(gene, p_val, avg_log2FC, abs_log2FC, pct_1, pct_2, p_val_adj);
    %    T = T(~isnan(p_val), :);
    %    T=sortrows(T,'p_val_adj','ascend');
    %    T=sortrows(T,'abs_logFC','descend');
    % [T] = pkg.e_sorttable(T);
    if nargout>1
        [Tup,Tdn]=pkg.e_processDETable(T);
    end
end

% Test for expression differences between two sets of cells
%
% Use the individual cell error models to test for differential expression between two groups of cells.
%
% @param models models determined by \code{\link{scde.error.models}}
% @param counts read count matrix
% @param prior gene expression prior as determined by \code{\link{scde.expression.prior}}
% @param groups a factor determining the two groups of cells being compared. The factor entries should correspond to the rows of the model matrix. The factor should have two levels. NAs are allowed (cells will be omitted from comparison).
% @param batch a factor (corresponding to rows of the model matrix) specifying batch assignment of each cell, to perform batch correction
% @param n.randomizations number of bootstrap randomizations to be performed
% @param n.cores number of cores to utilize
% @param batch.models (optional) separate models for the batch data (if generated using batch-specific group argument). Normally the same models are used.
% @param return.posteriors whether joint posterior matrices should be returned
% @param expectation M level corresponding to H0 hypothesis (usually 0, unless a deviation for a given gene is expected). Given on log2 scale.
% @param verbose integer verbose level (1 for verbose)
%
% @return \subsection{default}{
% a data frame with the following fields:
% \itemize{
% \item{lb, mle, ub} {lower bound, maximum likelihood estimate, and upper bound of the 95% confidence interval for the expression fold change on log2 scale.}
% \item{ce} { conservative estimate of expression-fold change (equals to the min(abs(c(lb, ub))), or 0 if the CI crosses the 0}
% \item{Z} { uncorrected Z-score of expression difference}
% \item{cZ} {expression difference Z-score corrected for multiple hypothesis testing using Benjamini-Hochberg procedure}
% }
%  If batch correction has been performed (\code{batch} has been supplied), analogous data frames are returned in slots \code{$batch.adjusted} for batch-corrected results, and \code{$batch.effect} for the differences explained by batch effects alone.
% }}
% \subsection{return.posteriors = TRUE}{
% A list is returned, with the default results data frame given in the \code{$results} slot.
% \code{difference.posterior} returns a matrix of estimated expression difference posteriors (rows - genes, columns correspond to different magnitudes of fold-change - log2 values are given in the column names)
% \code{joint.posteriors} a list of two joint posterior matrices (rows - genes, columns correspond to the expression levels, given by prior$x grid)
% }
%
% @examples
% data(es.mef.small)
% cd <- clean.counts(es.mef.small, min.lib.size=1000, min.reads = 1, min.detected = 1)
% sg <- factor(gsub("(MEF|ESC).*", "\\1", colnames(cd)), levels = c("ESC", "MEF"))
% names(sg) <- colnames(cd)
% \donttest{
% o.ifm <- scde.error.models(counts = cd, groups = sg, n.cores = 10, threshold.segmentation = TRUE)
% o.prior <- scde.expression.prior(models = o.ifm, counts = cd, length.out = 400, show.plot = FALSE)
% # make sure groups corresponds to the models (o.ifm)
% groups <- factor(gsub("(MEF|ESC).*", "\\1", rownames(o.ifm)), levels = c("ESC", "MEF"))
% names(groups) <- row.names(o.ifm)
% ediff <- scde.expression.difference(o.ifm, cd, o.prior, groups = groups, n.randomizations = 100, n.cores = n.cores, verbose = 1)
% }
%
% @export

%{
scde.expression.difference <- function(models, counts, prior, groups = NULL, batch = NULL, n.randomizations = 150, n.cores = 10, batch.models = models, return.posteriors = FALSE, expectation=0, verbose = 0) {
    if(!all(rownames(models) %in% colnames(counts))) {
        stop("ERROR: provided count data does not cover all of the cells specified in the model matrix")
    }

    ci <- match(rownames(models), colnames(counts))
    counts <- as.matrix(counts[, ci])

    # batch control
    if(correct.batch) {
        batch <- as.factor(batch)
        # check batch-group interactions
        bgti <- table(groups, batch)
        bgti.ft <- fisher.test(bgti)
        if(verbose) {
            cat("controlling for batch effects. interaction:\n")
            print(bgti)
        }
        #if(any(bgti == 0)) {
        #  cat("ERROR: cannot control for batch effect, as some batches are found only in one group:\n")
        #  print(bgti)
        #}
        if(bgti.ft$p.value < 1e-3) {
            cat("WARNING: strong interaction between groups and batches! Correction may be ineffective:\n")
            print(bgti.ft)
        }

        # calculate batch posterior
        if(verbose) {
            cat("calculating batch posteriors\n")
        }
        batch.jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
            scde.posteriors(models = batch.models, counts = counts, prior = prior, batch = batch, composition = table(batch[ii]), n.cores = n.cores, n.randomizations = n.randomizations, return.individual.posteriors = FALSE)
        })
        if(verbose) {
            cat("calculating batch differences\n")
        }
        batch.bdiffp <- calculate.ratio.posterior(batch.jpl[[1]], batch.jpl[[2]], prior, n.cores = n.cores)
        batch.bdiffp.rep <- quick.distribution.summary(batch.bdiffp)
    } else {
        if(verbose) {
            cat("comparing groups:\n")
            print(table(as.character(groups)))
        }
    }

    # fit joint posteriors for each group
    jpl <- tapply(seq_len(nrow(models)), groups, function(ii) {
        scde.posteriors(models = models[ii, , drop = FALSE], counts = counts[, ii, drop = FALSE], prior = prior, n.cores = n.cores, n.randomizations = n.randomizations)
    })
    if(verbose) {
        cat("calculating difference posterior\n")
    }

    # calculate difference posterior
    bdiffp <- calculate.ratio.posterior(jpl[[1]], jpl[[2]], prior, n.cores = n.cores)

    if(verbose) {
        cat("summarizing differences\n")
    }
    bdiffp.rep <- quick.distribution.summary(bdiffp,expectation=expectation)

    if(correct.batch) {
        if(verbose) {
            cat("adjusting for batch effects\n")
        }
        # adjust for batch effects
        a.bdiffp <- calculate.ratio.posterior(bdiffp, batch.bdiffp, prior = data.frame(x = as.numeric(colnames(bdiffp)), y = rep(1/ncol(bdiffp), ncol(bdiffp))), skip.prior.adjustment = TRUE, n.cores = n.cores)
        a.bdiffp.rep <- quick.distribution.summary(a.bdiffp,expectation=expectation)

        # return with batch correction info
        if(return.posteriors) {
            return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, batch.effect = batch.bdiffp.rep, difference.posterior = bdiffp, batch.adjusted.difference.posterior = a.bdiffp, joint.posteriors = jpl))
        } else {
            return(list(batch.adjusted = a.bdiffp.rep, results = bdiffp.rep, batch.effect = batch.bdiffp.rep))
        }
    } else {
        # no batch correction return
        if(return.posteriors) {
            return(list(results = bdiffp.rep, difference.posterior = bdiffp, joint.posteriors = jpl))
        } else {
            return(bdiffp.rep)
        }
    }
}

# calculates the likelihood of expression difference based on
# two posterior matrices (not adjusted for prior)
calculate.ratio.posterior <- function(pmat1, pmat2, prior, n.cores = 15, skip.prior.adjustment = FALSE) {
    n <- length(prior$x)
    if(!skip.prior.adjustment) {
        pmat1 <- t(t(pmat1)*prior$y)
        pmat2 <- t(t(pmat2)*prior$y)
    }

    chunk <- function(x, n) split(x, sort(rank(x) %% n))
    if(n.cores > 1) {
        x <- do.call(rbind, papply(chunk(1:nrow(pmat1), n.cores*5), function(ii) matSlideMult(pmat1[ii, , drop = FALSE], pmat2[ii, , drop = FALSE]), n.cores = n.cores))
    } else {
        x <- matSlideMult(pmat1, pmat2)
    }
    x <- x/rowSums(x)

    rv <- seq(prior$x[1]-prior$x[length(prior$x)], prior$x[length(prior$x)]-prior$x[1], length = length(prior$x)*2-1)
    colnames(x) <- as.character(rv)
    rownames(x) <- rownames(pmat1)
    return(x)
}

# given a set of pdfs (columns), calculate summary statistics (mle, 95% CI, Z-score deviations from 0)
# expectation - the M value representing H0 (usually 0), given on log2 scale
quick.distribution.summary <- function(s.bdiffp,expectation=0) {
    diffv <- as.numeric(colnames(s.bdiffp))
    dq <- t(apply(s.bdiffp, 1, function(p) {
        mle <- which.max(p)
        p <- cumsum(p)
        return(diffv[c(lb = max(c(1, which(p<0.025))), mle, min(c(length(p), which(p > (1-0.025)))))])
    }))/log10(2)
    colnames(dq) <- c("lb", "mle", "ub")
    cq <- rep(0, nrow(dq))
    cq[dq[, 1] > 0] <- dq[dq[, 1] > 0, 1]
    cq[dq[, 3]<0] <- dq[dq[, 3]<0, 3]
    z <- get.ratio.posterior.Z.score(s.bdiffp,expectation=expectation/log2(10))
    za <- sign(z)*qnorm(p.adjust(pnorm(abs(z), lower.tail = FALSE), method = "BH"), lower.tail = FALSE)
    data.frame(dq, "ce" = as.numeric(cq), "Z" = as.numeric(z), "cZ" = as.numeric(za))
}
https://raw.githubusercontent.com/hms-dbmi/scde/master/R/functions.R
%}

% fdr_bh() - Executes the Benjamini & Hochberg (1995) and the Benjamini &
%            Yekutieli (2001) procedure for controlling the false discovery
%            rate (FDR) of a family of hypothesis tests. FDR is the expected
%            proportion of rejected hypotheses that are mistakenly rejected
%            (i.e., the null hypothesis is actually true for those tests).
%            FDR is a somewhat less conservative/more powerful method for
%            correcting for multiple comparisons than procedures like Bonferroni
%            correction that provide strong control of the family-wise
%            error rate (i.e., the probability that one or more null
%            hypotheses are mistakenly rejected).
%
%            This function also returns the false coverage-statement rate
%            (FCR)-adjusted selected confidence interval coverage (i.e.,
%            the coverage needed to construct multiple comparison corrected
%            confidence intervals that correspond to the FDR-adjusted p-values).
%
%
% Usage:
%  >> [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(pvals,q,method,report);
%
% Required Input:
%   pvals - A vector or matrix (two dimensions or more) containing the
%           p-value of each individual test in a family of tests.
%
% Optional Inputs:
%   q       - The desired false discovery rate. {default: 0.05}
%   method  - ['pdep' or 'dep'] If 'pdep,' the original Bejnamini & Hochberg
%             FDR procedure is used, which is guaranteed to be accurate if
%             the individual tests are independent or positively dependent
%             (e.g., Gaussian variables that are positively correlated or
%             independent).  If 'dep,' the FDR procedure
%             described in Benjamini & Yekutieli (2001) that is guaranteed
%             to be accurate for any test dependency structure (e.g.,
%             Gaussian variables with any covariance matrix) is used. 'dep'
%             is always appropriate to use but is less powerful than 'pdep.'
%             {default: 'pdep'}
%   report  - ['yes' or 'no'] If 'yes', a brief summary of FDR results are
%             output to the MATLAB command line {default: 'no'}
%
%
% Outputs:
%   h       - A binary vector or matrix of the same size as the input "pvals."
%             If the ith element of h is 1, then the test that produced the
%             ith p-value in pvals is significant (i.e., the null hypothesis
%             of the test is rejected).
%   crit_p  - All uncorrected p-values less than or equal to crit_p are
%             significant (i.e., their null hypotheses are rejected).  If
%             no p-values are significant, crit_p=0.
%   adj_ci_cvrg - The FCR-adjusted BH- or BY-selected
%             confidence interval coverage. For any p-values that
%             are significant after FDR adjustment, this gives you the
%             proportion of coverage (e.g., 0.99) you should use when generating
%             confidence intervals for those parameters. In other words,
%             this allows you to correct your confidence intervals for
%             multiple comparisons. You can NOT obtain confidence intervals
%             for non-significant p-values. The adjusted confidence intervals
%             guarantee that the expected FCR is less than or equal to q
%             if using the appropriate FDR control algorithm for the
%             dependency structure of your data (Benjamini & Yekutieli, 2005).
%             FCR (i.e., false coverage-statement rate) is the proportion
%             of confidence intervals you construct
%             that miss the true value of the parameter. adj_ci=NaN if no
%             p-values are significant after adjustment.
%   adj_p   - All adjusted p-values less than or equal to q are significant
%             (i.e., their null hypotheses are rejected). Adjusted
%             p-values can be greater than 1.
%
%
% References:
%   Benjamini, Y. & Hochberg, Y. (1995) Controlling the false discovery
%     rate: A practical and powerful approach to multiple testing. Journal
%     of the Royal Statistical Society, Series B (Methodological). 57(1),
%     289-300.
%
%   Benjamini, Y. & Yekutieli, D. (2001) The control of the false discovery
%     rate in multiple testing under dependency. The Annals of Statistics.
%     29(4), 1165-1188.
%
%   Benjamini, Y., & Yekutieli, D. (2005). False discovery rate?adjusted
%     multiple confidence intervals for selected parameters. Journal of the
%     American Statistical Association, 100(469), 71?81. doi:10.1198/016214504000001907
%
%
% Example:
%  nullVars=randn(12,15);
%  [~, p_null]=ttest(nullVars); %15 tests where the null hypothesis
%  %is true
%  effectVars=randn(12,5)+1;
%  [~, p_effect]=ttest(effectVars); %5 tests where the null
%  %hypothesis is false
%  [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh([p_null p_effect],.05,'pdep','yes');
%  data=[nullVars effectVars];
%  fcr_adj_cis=NaN*zeros(2,20); %initialize confidence interval bounds to NaN
%  if ~isnan(adj_ci_cvrg),
%     sigIds=find(h);
%     fcr_adj_cis(:,sigIds)=tCIs(data(:,sigIds),adj_ci_cvrg); % tCIs.m is available on the
%     %Mathworks File Exchagne
%  end
%
%
% For a review of false discovery rate control and other contemporary
% techniques for correcting for multiple comparisons see:
%
%   Groppe, D.M., Urbach, T.P., & Kutas, M. (2011) Mass univariate analysis
% of event-related brain potentials/fields I: A critical tutorial review.
% Psychophysiology, 48(12) pp. 1711-1725, DOI: 10.1111/j.1469-8986.2011.01273.x
% http://www.cogsci.ucsd.edu/~dgroppe/PUBLICATIONS/mass_uni_preprint1.pdf
%
%
% For a review of FCR-adjusted confidence intervals (CIs) and other techniques
% for adjusting CIs for multiple comparisons see:
%
%   Groppe, D.M. (in press) Combating the scientific decline effect with
% confidence (intervals). Psychophysiology.
% http://biorxiv.org/content/biorxiv/early/2015/12/10/034074.full.pdf
%
%
% Author:
% David M. Groppe
% Kutaslab
% Dept. of Cognitive Science
% University of California, San Diego
% March 24, 2010

%%%%%%%%%%%%%%%% REVISION LOG %%%%%%%%%%%%%%%%%
%
% 5/7/2010-Added FDR adjusted p-values
% 5/14/2013- D.H.J. Poot, Erasmus MC, improved run-time complexity
% 10/2015- Now returns FCR adjusted confidence intervals

function [h, crit_p, adj_ci_cvrg, adj_p] = fdr_bh(pvals, q, method, report)

    if nargin < 1
        error('You need to provide a vector or matrix of p-values.');
    else
        if ~isempty(find(pvals < 0, 1))
            error('Some p-values are less than 0.');
        elseif ~isempty(find(pvals > 1, 1))
            error('Some p-values are greater than 1.');
        end
    end

    if nargin < 2
        q = .05;
    end

    if nargin < 3
        method = 'pdep';
    end

    if nargin < 4
        report = 'no';
    end

    s = size(pvals);
    if (length(s) > 2) || s(1) > 1
        [p_sorted, sort_ids] = sort(reshape(pvals, 1, prod(s)));
    else
        % p-values are already a row vector
        [p_sorted, sort_ids] = sort(pvals);
    end
    [~, unsort_ids] = sort(sort_ids); % indexes to return p_sorted to pvals order
    m = length(p_sorted); % number of tests

    if strcmpi(method, 'pdep')
        % BH procedure for independence or positive dependence
        thresh = (1:m) * q / m;
        wtd_p = m * p_sorted ./ (1:m);

    elseif strcmpi(method, 'dep')
        % BH procedure for any dependency structure
        denom = m * sum(1 ./ (1:m));
        thresh = (1:m) * q / denom;
        wtd_p = denom * p_sorted ./ (1:m);
        % it can produce adjusted p-values greater than 1!
        % compute adjusted p-values
    else
        error('Argument ''method'' needs to be ''pdep'' or ''dep''.');
    end

    if nargout > 3
        % compute adjusted p-values; This can be a bit computationally intensive
        adj_p = zeros(1, m) * NaN;
        [wtd_p_sorted, wtd_p_sindex] = sort(wtd_p);
        nextfill = 1;
        for k = 1:m
            if wtd_p_sindex(k) >= nextfill
                adj_p(nextfill:wtd_p_sindex(k)) = wtd_p_sorted(k);
                nextfill = wtd_p_sindex(k) + 1;
                if nextfill > m
                    break
                end
            end
        end
        adj_p = reshape(adj_p(unsort_ids), s);
    end

    rej = p_sorted <= thresh;
    max_id = find(rej, 1, 'last'); % find greatest significant pvalue
    if isempty(max_id)
        crit_p = 0;
        h = pvals * 0;
        adj_ci_cvrg = NaN;
    else
        crit_p = p_sorted(max_id);
        h = pvals <= crit_p;
        adj_ci_cvrg = 1 - thresh(max_id);
    end

    if strcmpi(report, 'yes')
        n_sig = sum(p_sorted <= crit_p);
        if n_sig == 1
            fprintf('Out of %d tests, %d is significant using a false discovery rate of %f.\n', m, n_sig, q);
        else
            fprintf('Out of %d tests, %d are significant using a false discovery rate of %f.\n', m, n_sig, q);
        end
        if strcmpi(method, 'pdep')
            fprintf('FDR/FCR procedure used is guaranteed valid for independent or positively dependent tests.\n');
        else
            fprintf('FDR/FCR procedure used is guaranteed valid for independent or dependent tests.\n');
        end
    end
end
