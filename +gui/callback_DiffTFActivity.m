function callback_DiffTFActivity(src,~)
    
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    [thisc,clable]=gui.i_select1class(sce,false);
    species=gui.i_selectspecies(2);

    fw=gui.gui_waitbar;
    [cs,tflist,~,numtargetgenes]=sc_tfactivity(sce.X,sce.g,[],species);    
    T=pkg.e_grptest(cs,thisc);
    T=[table(tflist), T, table(numtargetgenes)];
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

    [~,cL]=grp2idx(thisc);
    cL=strrep(cL,'_','\_');
    thisc=strrep(thisc,'_','\_');
    
    answer=questdlg('Violin plot for top 10 TFs with most variable activity levels between groups?');

    switch answer
        case 'Yes'
%             a=inputdlg('Number of top TFs:','',1,{'10'});
%             if isempty(a), return; end
%             a=str2double(a{1});
            [numfig]=gui.i_inputnumg(length(T.tflist));

            if isempty(numfig) || isnan(numfig), return; end
            if isnumeric(numfig) && numfig>0 && numfig<=length(T.tflist)
                F=cell(numfig,1);
                for k=1:numfig        
                    f = figure('visible','off');
                    pkg.i_violinplot(cs(k,:),thisc,true,cL);
                    title(T.tflist(k));
                    ylabel('TF activity')
                    P = get(f,'Position');
                    set(f,'Position',[P(1)-20*k P(2)-20*k P(3) P(4)]);
                    set(f,'visible','on');
                    f.Position(3) = f.Position(3) * 2.2;
                    drawnow;
                    F{k}=f;
                end
                gui.i_export2pptx(F,string(T.tflist(1:numfig)));
            end
    end
end
