% [c, W] = run_soptsc(X, 3);
% 
% cluster_label = c;
% No_cluster = length(unique(c));
% root_cluster = 0;
% root_cell = 0;
% reverse = 1;
% 
% W1 = W ./ (ones(1, size(W, 1)) * W * ones(size(W, 1), 1));
% InitY = pca1(W1, 2);
% latent = tsne(W1, 'Standardize', true, 'Perplexity', 35, 'NumDimensions', 2, 'InitialY', InitY);
% 
% [Lineage, Ptime, Cell_dist] = Lineage_Ptime(W, No_cluster, cluster_label, root_cluster, root_cell, latent, reverse);
% 
% plot_lineage(Lineage, No_cluster, cluster_label, Cell_dist)
