function callback_DEGene2Groups(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    

    listitems={'Cluster ID','Batch ID',...
        'Cell Type','Cell Cycle Phase'};
    [indx2,tf2] = listdlg('PromptString',{'Select statistics',...
    '',''},'SelectionMode','single','ListString',listitems);
    if tf2==1        
        switch indx2
            case 1 % cluster id
                thisc=sce.c_cluster_id;
            case 2 % batch id
                thisc=sce.c_batch_id;
            case 3 % cell type
                thisc=sce.c_cell_type_tx;                
            case 4 % cell cycle
                thisc=sce.c_cell_cycle_tx;
        end
    else
        return;
    end 
    if isempty(thisc)
        errordlg('Undefined');
        return;
    end
    if numel(unique(thisc))==1
        warndlg("Cannot compare with an unique group");
        return;
    end
    
    [ci,cLi]=grp2idx(thisc);
    [indxx,tfx] = listdlg('PromptString',{'Select two groups',...
    '',''},'SelectionMode','multiple','ListString',string(cLi),...
    'InitialValue',[1 2]);
    if tfx==1
        if numel(indxx)~=2
            errordlg('Please select 2 groups');
            return;
        end
        i1=ismember(ci,indxx(1));
        i2=ismember(ci,indxx(2));
    else
        return;
    end

    

    answer = questdlg('Which method?','Select Method','Wilcoxon rank-sum test','MAST','Wilcoxon rank-sum test');
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
    labels = {'Save DE results T to variable named:'}; 
    vars = {'T'}; values = {T};
    msgfig=export2wsdlg(labels,vars,values);
    uiwait(msgfig);
    answer = questdlg('Violin plots?');
    if strcmp(answer,'Yes')
        figure;
        for k=1:16
            subplot(4,4,k)
            i=sce.g==T.gene(k);
            pkg.i_violinplot(log2(1+sce.X(i,:)),...
                sce.c_batch_id);
            title(T.gene(k));
            ylabel('log2(UMI+1)')
        end
    end
    

            
end