function plot_classic(ES,ESmax,data_ranks,RI_GS,name,Phit,Pmiss)
% Plot GSEA results for each gene set

if isnumeric(name)
    name = num2str(name);   
else
    name = char(name);
    name(name=='_')=' '; % for names with _ removing _
end

f=figure('Visible','off');
subplot(3,1,1); hold on; box on
plot(1:length(ES),ES,'b','LineWidth',1.5)
plot(repmat(find(ES==ESmax),[1 2]),[0,ESmax],'--r',1:length(ES),zeros(1,length(ES)),'k');
xlim([1,length(ES)]); ylabel('Enrichment Score')
title(['GeneSet ID: ',name])
subplot(3,1,2); hold on; box on
plot(cumsum(Pmiss),'b-','LineWidth',2); plot(cumsum(Phit),'-r','LineWidth',2);   
ylabel('Cumulative distribution');
legend({'Pmiss','Phit'},'Location','southeast')
ylim([0 1]); xlim([1 length(ES)])
subplot(3,1,3); hold on; box on
bar(data_ranks)
plot([RI_GS; RI_GS],[zeros(1,length(RI_GS)); repmat(max(data_ranks),[1 length(RI_GS)])],'-r')
xlim([1 length(ES)]); ylabel('Rank distribution'); xlabel('Gene order')
% saveas(gcf,['GSEA_plots/GS_',name,'.png'],'png');
saveas(f,['GSEA_plots/GS_',name,'.png'],'png');

