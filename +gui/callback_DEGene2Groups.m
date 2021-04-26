function callback_DEGene2Groups(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    

    [i1,i2]=gui.i_select2grps(sce);
    if length(i1)==1 || length(i2)==1, return; end

    answer = questdlg('Which method?',...
        'Select Method','Wilcoxon rank-sum test',...
        'MAST','Wilcoxon rank-sum test');
    
    if strcmpi(answer,'Wilcoxon rank-sum test')
        methodtag="ranksum";
    elseif strcmpi(answer,'MAST')
        methodtag="mast";
    else
        return;
    end
    fw=gui.gui_waitbar;
    switch methodtag
        case 'ranksum'
            T=sc_deg(sce.X(:,i1),...
                    sce.X(:,i2),sce.g);
        case 'mast'
            T=run.MAST(sce.X(:,i1),...
                    sce.X(:,i2),sce.g);
    end
    gui.gui_waitbar(fw);
    isok=(abs(T.pct_2-T.pct_1)>0.05 | abs(T.avg_logFC)>1.0)&T.p_val_adj<0.01;
    Tup=T(T.avg_logFC>0 & isok,:);
    Tdn=T(T.avg_logFC<0 & isok,:);
    
Tup = sortrows(Tup,'abs_logFC','descend');
Tup = sortrows(Tup,'p_val_adj','ascend');

Tdn = sortrows(Tdn,'abs_logFC','descend');
Tdn = sortrows(Tdn,'p_val_adj','ascend');

    labels = {'Save DE results (all genes) to variable named:',...
        'Save DE results (selected up-regulated) to variable named:',...
        'Save DE results (selected down-regulated) to variable named:'}; 
    vars = {'T','Tup','Tdn'}; values = {T,Tup,Tdn};
    msgfig=export2wsdlg(labels,vars,values);
    uiwait(msgfig);
    
    answer = questdlg('Run Enrichr with top 200 up-regulated DE genes?');
    if strcmp(answer,'Yes')
        run.Enrichr(Tup.gene(1:min(numel(Tup.gene),200)))
    end
    pause(3);
    answer = questdlg('Run Enrichr with top 200 down-regulated DE genes?');
    if strcmp(answer,'Yes')
        run.Enrichr(Tdn.gene(1:min(numel(Tdn.gene),200)))
    end
    pause(3);
    
    %Tf=run.fgsea(T.gene);
    %pkg.e_fgseanet(Tf);
    
%     answer = questdlg('Violin plots (top 16 DE genes)?');
%     if strcmp(answer,'Yes')
%         figure;
%         for k=1:16
%             subplot(4,4,k)
%             i=sce.g==T.gene(k);
%             pkg.i_violinplot(log2(1+sce.X(i,:)),...
%                 sce.c_batch_id);
%             title(T.gene(k));
%             ylabel('log2(UMI+1)')
%         end
%     end

    answer = questdlg('Volcano plot?');
    if strcmp(answer,'Yes')
        figure;
        gui.gui_volcanoplot(T,isok);
    end
end