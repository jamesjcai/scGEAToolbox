function callback_CompareGeneBtwCls(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    

    gsorted=sort(sce.g);
    [indx,tf] = listdlg('PromptString',{'Select a gene',...
    '',''},'SelectionMode','single','ListString',gsorted);
    if tf==1
        idx=sce.g==gsorted(indx);
%{        
        answer = questdlg('Which category?','',...
            'Cluster ID','Cell Type','Cell Cycle','Cluster ID');
        if strcmp(answer,'Cluster ID')
            thisc=sce.c_cluster_id;
        elseif strcmp(answer,'Cell Type')
            thisc=sce.c_cell_type_tx;
        elseif strcmp(answer,'Cell Cycle')
            thisc=sce.c_cell_cycle_tx;
        else
            return;
        end
%}
        
    listitems={'Cluster ID','Batch ID',...
        'Cell Type','Cell Cycle Phase'};
    [indx2,tf2] = listdlg('PromptString',{'Select statistics',...
    '',''},'SelectionMode','single','ListString',listitems);
    if tf2==1        
        switch indx2
            case 1
                thisc=sce.c_cluster_id;
            case 2
                thisc=sce.c_batch_id;
            case 3
                thisc=sce.c_cell_type_tx;                
            case 4
                thisc=sce.c_cell_cycle_tx;
        end
    else
        return;
    end 
        if isempty(thisc)
            errordlg('Undefined');
        else
            [Xt]=gui.gui_transformx(X);      
            f = figure('visible','off');
            y=Xt(idx,:);
            subplot(1,2,1)
            pkg.violinplot(y,thisc);
            ylabel(sce.g(idx));
            subplot(1,2,2); 
            hold on
            [a,b]=grp2idx(thisc);
            for kkx=1:max(a)
                cdfplot(y(a==kkx));
            end
            legend(b);
            f.Position(3)=round(f.Position(3)*1.5);           
            f.Position(4)=round(f.Position(4)*0.7);
            movegui(f,'center');
            set(f,'visible','on');
        end
    end
end