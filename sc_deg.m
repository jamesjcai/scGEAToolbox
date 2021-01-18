function [T]=sc_deg(X,Y,genelist,methodid)
% https://satijalab.org/seurat/v3.1/de_vignette.html
% p_val : p_val (unadjusted)
% avg_logFC : log fold-chage of the average expression between the two groups. Positive values indicate that the feature is more highly expressed in the first group.
% pct.1 : The percentage of cells where the feature is detected in the first group
% pct.2 : The percentage of cells where the feature is detected in the second group
% p_val_adj : Adjusted p-value, based on bonferroni correction using all features in the dataset.
%
% SEE ALSO: [T]=run_mast(X,Y,genelist);

if nargin<2, error("USAGE: sc_deg(X,Y)\n"); end
if nargin<3, genelist=string(num2cell(1:size(X,1)))'; end
if nargin<4, methodid=1; end

ng=size(X,1);
assert(isequal(ng,size(Y,1)));

p_val=ones(ng,1);
avg_logFC=ones(ng,1);
pct_1=ones(ng,1);
pct_2=ones(ng,1);

for k=1:ng
    x=X(k,:);
    y=Y(k,:);    
    switch methodid
        case 1
            % “wilcox” : Wilcoxon rank sum test (default)
            % also called Mann–Whitney U test
            p_val(k)=ranksum(x,y);
        case 2
            [~,p]=ttest2(x,y);
            p_val(k)=p;
        otherwise
            p_val(k)=ranksum(x,y);
    end
    avg_logFC(k)=log2(mean(x)./mean(y));
    pct_1(k)=sum(x>0)./length(x);
    pct_2(k)=sum(y>0)./length(y);
end

    p_val_adj = mafdr(p_val,'BHFDR',true);

%     sortid=(1:length(genelist))';
    if size(genelist,2)>1 
        gene=genelist';
    else
        gene=genelist;
    end
    T=table(gene,p_val,avg_logFC,pct_1,pct_2,p_val_adj);
    T=T(~isnan(p_val),:);
    T=sortrows(T,'p_val_adj','ascend');


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
