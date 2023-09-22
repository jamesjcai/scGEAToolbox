function plot_marker(data, gene_set, allgenes, latent)
%This function plot gene expression along cell subpopulations
%   Input:
%  --   gene_set: genes to be plot specified by users
%  --   all_genes: all genes from the single cell data set
%  --   data: a m*n single cell data matrix with m rows(genes) and n columns(cells)
%  --   latent: low dimension space induced from cell-cell transition matrix
%  --   folder: folder name where the results will be save to.
%
%   Output:
%           figures showing gene expression among cell subpopulations
%           on 2D-projection.
%

gene_set_no = length(gene_set);
allgenes_no = length(allgenes);


gene_set_idx = zeros(1, gene_set_no);

for i = 1:gene_set_no
    for j = 1:allgenes_no
        if strcmpi(gene_set{i}, allgenes{j})
            gene_set_idx(i) = j;
        end
    end
end
%gene_set_idx

MM = data(gene_set_idx, :);

%% Marker genes expression on each subpopulation

% load newmap.mat;
% mymap = newmap;

for ik = 1:gene_set_no
    figure;
    colormap(jet);
    scatter(latent(:, 1), latent(:, 2), 20, MM(ik, :), 'filled', 'MarkerEdgeAlpha', 0.8, 'MarkerFaceAlpha', 0.8);
    title(gene_set{ik});
    box off;
    %     set(gca,'LineWidth',1.5);
    set(gca, 'xtick', []);
    set(gca, 'ytick', []);
    set(gca, 'FontName', 'Arial');
    set(gca, 'FontSize', 12);

    ax = gca;
    cb = colorbar;
    % ax = gca;
    axpos = ax.Position;
    cpos = cb.Position;
    cpos(3) = 0.5 * cpos(3);
    cb.Position = cpos;
    ax.Position = axpos;
    % lim = caxis
    % cb.Limits = lim;
    aa = cell(1, 2);
    aa{1} = 'low';
    aa{2} = 'high';
    cb.TickLabels{1} = aa{1};
    cb.TickLabels{end} = aa{2};

    for ii = 2:length(cb.TickLabels) - 1
        cb.TickLabels{ii} = [];
    end
    % print([resfolder '\mk_ldm_' gene_set{ik}],'-dpdf','-r300');
end
end
