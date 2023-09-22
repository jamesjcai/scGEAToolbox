function i_qcviolin(X, genelist)

i = startsWith(genelist, 'mt-', 'IgnoreCase', true);
nftr = full(sum(X > 0, 1));
lbsz = full(sum(X, 1));
lbsz_mt = full(sum(X(i, :), 1));
cj = 100 * (lbsz_mt ./ lbsz);


figure;
subplot(1, 3, 1)
pkg.violinplot(nftr, [], 'showdata', false);
title(sprintf('nFeature\\_RNA\n(# of genes)'));
box on;

subplot(1, 3, 2)
pkg.violinplot(lbsz, [], 'showdata', false);
title(sprintf('nCount\\_RNA\n(# of reads)'));
box on;

subplot(1, 3, 3)
pkg.violinplot(cj, [], 'showdata', false);
title(sprintf('percent.mt\n(mitochondrial content)'));
box on;
