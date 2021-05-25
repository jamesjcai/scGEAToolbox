function opts = default_GSEA_opts()
%Default parameters for GSEA method.

opts.GS_name = 'KEGG_GS';   %GeneSet database name
opts.GS_filt = [15,500];    %minimum and maximum number of genes in GS
opts.sort_type = 'descend'; %type of ranks sorting
opts.p = 1;                 %0 - Kolmogorov-Smirnov statistic, 1-weighting genes by ranking metric
opts.rank_type = 'S2N';     %ranking method {'S2N','ttest','ratio','diff','log2_ratio','SoR','BWS','ReliefF','WAD','FCROS','MWT'}
opts.abs = true;            %if use absolute value of ranking statistic
opts.tied_rank = false;     %if create tied ranks of ranking statistic
opts.perm_type = 'pheno';   %type of permutation {'entrez','pheno'}
opts.perm_nb = 1000;        %number of permutations
opts.perm_ext = false;      %if extend number of permutation till one significant result occur  (may be time consuming)
opts.GS_sort_type = 'none'; %options for sorting results {'ES_pval','NES_qval','NES_FWER','none'}
opts.show = true;           %if plot results
opts.save = false;          %if save results
