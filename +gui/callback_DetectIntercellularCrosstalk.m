function callback_DetectIntercellularCrosstalk(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
    if isempty(sce.c_cell_type_tx) || numel(unique(sce.c_cell_type_tx))<2
        if ~isempty(sce.c_cluster_id) && numel(unique(sce.c_cluster_id))>1
            answer = questdlg(sprintf('Cell type (C_CELL_TYPE_TX) is undefined.\nWould you like to use cluster id (C_CLUSTER_ID) to define cell groups?'));
            switch answer
                case 'Yes'
                    sce.c_cell_type_tx=strcat('Goup', string(sce.c_cluster_id));
                otherwise
                    return;
            end
        else
            warndlg('Cell type is undefined (SCE.C_CELL_TYPE_TX is empty)');
            return;
        end
    end
    
    
    [c,cL]=grp2idx(sce.c_cell_type_tx);
    [idx]=gui.i_selmultidlg(cL);
    if isempty(idx), return; end
    if numel(idx)<2
        warndlg('Need at least 2 cell types');
        return;
    end
    selected=ismember(c,idx);
    fw=gui.gui_waitbar;
    sce=sce.selectcells(selected);
    [OUT]=run.talklr(sce);
    gui.gui_waitbar(fw);    
    
    n=length(OUT.ligandok);
    if n==0
        warndlg('Not detected.');
    end
    
    labels = {'Save OUT to variable named:'};
    vars = {'OUT'};
    values = {OUT};
    [f,ft]=export2wsdlg(labels, vars, values);
    waitfor(f);
    if ft
        disp('Run gui.i_crosstalkgraph(OUT,k) to plot crosstalk graph for ligand-receptor pair k.')
    end



    listitems=cell(n,1);
    for k=1:n
        listitems{k}=sprintf('%s -> %s (KL = %.2f)',...
            OUT.ligandok(k), OUT.receptorok(k),...
            OUT.KL(k));
    end
    [indx2,tf2] = listdlg('PromptString',...
        {'Select ligand-receptor pairs to plot','',''},...
         'SelectionMode','multiple','ListString',listitems,...
         'ListSize',[210,300]);
    if tf2==1
        for kk=1:length(indx2)
                [y,idx]=ismember(upper(OUT.ligandok(kk)),upper(sce.g));
                if y
                    figure;
                    sc_scattermarker(sce.X,sce.g,sce.s,sce.g(idx),5);
                end                
                [y,idx]=ismember(upper(OUT.receptorok(kk)),upper(sce.g));
                if y
                    figure;
                    sc_scattermarker(sce.X,sce.g,sce.s,sce.g(idx),5);
                end
                gui.i_crosstalkgraph(OUT,kk);
        end
    end    
end
