function latent = plot_cluster(W,cluster_label,No_cluster,method,~)

% plot cluster on 2-dimensional space
% input
% - W: similarity matrix
% - eigenvalues: eigenvalues of graph Laplacian of the truncated consensus
%   matrix
% - method: dimension reduction method, 'tsne' or 'pca' (e.g., method = 'tsne')
%
% Output
% - latent: 2-dimensional projection

W1 = W./(ones(1,size(W,1))*W*ones(size(W,1),1));
if strcmp(method,'pca')
    dvis = pca1(W1,2);
elseif strcmp(method,'tsne')
        InitY = pca1(W1,2);
        dvis = tsne(W1,'Standardize',true,'Perplexity',35,'NumDimensions',2,'InitialY',InitY);
end


%% Subpopulations visualization
% colormap jet;
cmap1 = jet;
mymap1 = cmap1(1:end,:);
ncolor = size(mymap1,1);

mycolor = mymap1(1:round(ncolor./No_cluster):ncolor,:);

figure;
for ik = 1:No_cluster
    scatter(dvis(find(cluster_label==ik),1),dvis(find(cluster_label==ik),2),40,mycolor(ik,:),'filled','MarkerEdgeAlpha',0.6,'MarkerFaceAlpha',0.6);
    hold on;
end
box off;
%     set(gca,'LineWidth',1.5);
set(gca,'xtick',[]);
set(gca,'ytick',[]);

lgd = cell(1,No_cluster);
for i = 1:No_cluster
    if i<10
        vv = 'CC';
        vv(2:2) = num2str(i);
        lgd{i} = vv;
    else
        vv = 'CCC';
        vv(2:3) = num2str(i);
        lgd{i} = vv;
    end
end
legend(lgd,'FontSize',12,'Location','best');%,'Orientation','horizontal');
set(gca,'FontName','Arial');
set(gca,'FontSize',12);


ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 4 3];
%
%     fig.Units = 'Inches';
%     fig.Position = [0 0 4 3];


%print([folder '\Cluster_' method],'-dpdf','-r300','-fillpage'); %'-dpdf',
latent = dvis;
