function callback_DiffTFActivity(src,~)
    
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc,clable]=gui.i_select1class(sce,false);
    species=gui.i_selectspecies;

    fw=gui.gui_waitbar;
    [cs,tflist]=sc_tfactivity(sce.X,sce.g,[],species);    
    T=pkg.e_grptest(cs,thisc);
    T=[table(tflist), T];
    gui.gui_waitbar(fw);

    
    [yis]=ismember(upper(tflist),upper(sce.g));
    T2=T(yis,:);
    cs2=cs(yis,:);

    outfile=sprintf('DiffTFActivity_%s',matlab.lang.makeValidName(clable));
    
    answer=questdlg(sprintf('%d out of %d TFs expressed in cells. Keep only %d expressed TFs?',...
        size(T2,1),size(T,1),size(T2,1)));
    switch answer
        case 'Yes'
            T=T2;
            cs=cs2;
        case 'No'          
        otherwise
            return;
    end

    try
        if length(unique(thisc))>2
            [T,idx] = sortrows(T,"p_anova","ascend");
            cs=cs(idx,:);
            [T,idx] = sortrows(T,"p_kruskalwallis","ascend");            
            cs=cs(idx,:);
        else
            [T,idx] = sortrows(T,"p_ttest","ascend");
            cs=cs(idx,:);
            [T,idx] = sortrows(T,"p_wilcoxon","ascend");
            cs=cs(idx,:);
        end
    catch

    end
    gui.i_exporttable(T,true,'T',outfile);

    for k=1:10
        [~,cL]=grp2idx(thisc);
        cL=strrep(cL,'_','\_');
        thisc=strrep(thisc,'_','\_');
        figure;
        pkg.i_violinplot(cs(k,:),thisc,true,cL);
        title(T.tflist(k));
    end
end
