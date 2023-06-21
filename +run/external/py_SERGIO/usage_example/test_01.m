%load('opsahl_11.mat','x','comm_real')
%coords = coalescent_embedding(x, 'RA1', 'ISO', 'original', 3);
%[coords(:,1),coords(:,2),coords(:,3)] = sph2cart(coords(:,1),coords(:,2),coords(:,3));
%coords=coords./vecnorm(coords,2,2);

load('PSO_network.mat', 'x', 'coords_orig')
coords_emb = coalescent_embedding(x, 'RA1', 'LE', 'original', 2);

coords=coords_emb;
[coords(:,1),coords(:,2)] = pol2cart(coords_emb(:,1),coords_emb(:,2));
coords=coords./vecnorm(coords,2,2);
figure; scatter(coords(:,1),coords(:,2),[],coords_orig(:,1));

%%

% evaluate comparison of embedding with original coordinates
[HD_corr, C_score] = compare_embedding(coords_orig, coords_emb); %#ok<NODEF>

% plot comparison of embedding with original coordinates
figure('color','white')
coords_emb_aligned = coords_emb;
coords_emb_aligned(:,1) = angular_alignment(coords_orig(:,1), coords_emb(:,1), 360);

subplot(1,3,1)
plot_embedding(x, coords_orig, 'labels', coords_orig(:,1));
text(0.50, 1.1, 'Original coordinates', 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 20)
subplot(1,3,2)
plot_embedding(x, coords_emb_aligned, 'labels', coords_orig(:,1));
text(0.50, 1.1, 'Embedding coordinates (RA1-LE-EA)', 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 20)
text(0.50, -0.1, ['HD-correlation = ' num2str(HD_corr,2)], 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 16)
subplot(1,3,3)
scatter(coords_orig(:,1), coords_emb_aligned(:,1), 'filled', 'b');
xlabel('Original angular coordinates', 'FontSize', 14); ylabel('Embedding angular coordinates', 'FontSize', 14);
xlim([0,2*pi]); ylim([0,2*pi]); axis square; box on;
text(0.50, 1.1, ['C-score = ' num2str(C_score,2)], 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 16)
   