% library(apcluster)
% lk <- linKernel(MA, normalize = TRUE)
% outCluster <- apcluster(s = lk[1:100,1:100])
% save(outCluster, file = 'outCluster.RData')

[~, i] = sort(vecnorm(aln0-aln1, 2, 2), 'descend');
gx = upper(unique(genelist(i), 'stable'));

