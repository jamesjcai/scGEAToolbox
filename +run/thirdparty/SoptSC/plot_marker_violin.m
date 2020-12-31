function plot_marker_violin(data,allgenes,Marker,cluster_label,No_cluster)

pw1=fileparts(which(mfilename));
pth=fullfile(pw1,'vinlinplot');
addpath(pth);

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



cmap1 = jet;
mymap1 = cmap1(1:end,:);
ncolor = size(mymap1,1);
mycolor = mymap1(1:round(ncolor./No_cluster):ncolor,:);


cluster_label1 = [];
cell_order = [];
for i = 1:No_cluster
    cell_order = [cell_order; find(cluster_label==i)];
    cluster_label1 = [cluster_label1; i*ones(length(find(cluster_label==i)),1)];
end


for i = 1:length(Marker)
[~,ia,~] = intersect(allgenes,Marker{i},'stable');
figure
h = violinplot(data(ia,cell_order),cluster_label1,1:No_cluster);
    for j = 1:No_cluster
        h(j).ViolinColor = mycolor(j,:);
    end
set(gca,'Xtick',1:No_cluster)
set(gca,'Xticklabel',lgd,'FontSize',12);
title(Marker{i});
% ax = gca;
% outerpos = ax.OuterPosition;
% ti = ax.TightInset; 
% left = outerpos(1) + ti(1);
% bottom = outerpos(2) + ti(2);
% ax_width = outerpos(3) - ti(1) - ti(3);
% ax_height = outerpos(4) - ti(2) - ti(4);
% ax.Position = [left bottom ax_width ax_height];
% 
% fig = gcf;
% fig.PaperPositionMode = 'auto';
% fig_pos = fig.PaperPosition;
% fig.PaperSize = [fig_pos(3) fig_pos(4)];
% 
% ax.XTickMode = 'manual';
% ax.YTickMode = 'manual';
% ax.ZTickMode = 'manual';
% ax.XLimMode = 'manual';
% ax.YLimMode = 'manual';
% ax.ZLimMode = 'manual';

box off;
grid on;


ax = gca;
fig = gcf;
fig.PaperPositionMode = 'auto';
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

ax.XTickMode = 'manual';
ax.YTickMode = 'manual';
ax.ZTickMode = 'manual';
ax.XLimMode = 'manual';
ax.YLimMode = 'manual';
ax.ZLimMode = 'manual';

% print([folder '\mk_violin_' Marker{i}],'-dpdf','-r300'); %'-dpdf',
end