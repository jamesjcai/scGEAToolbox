function RUN_EXAMPLE(example_num)

% Authors:
% - main code: Alessandro Muscoloni, 2017-09-21

% Released under MIT License
% Copyright (c) 2017 A. Muscoloni, J. M. Thomas, C. V. Cannistraci

% Usage:
% RUN_EXAMPLE(1)
% RUN_EXAMPLE(2)

% Example 1:
% - 2D embedding of PSO network using RA1-LE-EA
% - evaluate comparison of embedding with original coordinates
% - plot comparison of embedding with original coordinates
% Example 2:
% - 3D embedding of opsahl-11 network using RA1-ISO
% - evaluation of 3D greedy routing
% - 3D plot colored by ground-truth communities

narginchk(1,1)
if example_num == 1
    
    %%% Example 1 %%%
    
    % 2D embedding of PSO network using RA1-LE-EA
    load('PSO_network.mat', 'x', 'coords_orig')
    display('2D embedding of PSO network using RA1-LE-EA')
    coords_emb = coalescent_embedding(x, 'RA1', 'LE', 'EA', 2);
    
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
        
elseif example_num == 2
    
    %%% Example 2 %%%
    
    % 3D embedding of opsahl-11 network using RA1-ISO
    load('opsahl_11.mat','x','comm_real')
    display('3D embedding of opsahl-11 network using RA1-ISO')
    coords = coalescent_embedding(x, 'RA1', 'ISO', 'original', 3);

    
    % evaluation of 3D greedy routing
    display('Greedy routing...')
    [GR_score, succ_paths, avg_length] = greedy_routing(x, coords);
    
    % 3D plot colored by ground-truth communities
    figure('color','white')
    plot_embedding(x, coords, 'labels', comm_real);
    text(0.50, 1.05, '3D embedding of opsahl-11 network using RA1-ISO', 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 20)
    text(0.50, -0.05, ['Greedy routing: GR-score = ' num2str(GR_score,2) '; succ-paths = ' num2str(succ_paths,2) '; avg-length = ' num2str(avg_length,2)], 'horizontalAlignment', 'center', 'Units', 'normalized', 'FontSize', 20)
        
else
	error('Usage: RUN_EXAMPLE(1) or RUN_EXAMPLE(2).')
end
