function callback_ScatterStemPlot(src, ~)

    if isa(src,"SingleCellExperiment")
        sce = src;
        FigureHandle = [];
    else
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end
    
    
    answer = gui.myQuestdlg(FigureHandle, "Scatter-stem plot for gene expression or cell state variables?","", ...
        {'Gene Expression', 'Cell State'}, 'Gene Expression');
    switch answer
        case 'Gene Expression'
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist)
                gui.myHelpdlg(FigureHandle, 'No gene selected.', '');
                return;
            end

            [Xt] = gui.i_transformx(sce.X, [], [], FigureHandle);
            if isempty(Xt), return; end

            % answer = gui.myQuestdlg(FigureHandle, 'Plot all in the same figure?','');
            % if strcmp(answer, 'Yes')                
            %     fw = gui.myWaitbar(FigureHandle);                
            %     gui.i_violinmatrix(full(Xt), sce.g, c, cL, glist, ...
            %             false, '', FigureHandle);
            % 
            %     gui.myWaitbar(FigureHandle, fw);
            % 
            %     return;
            % end

            n = length(glist);
            thisyv = cell(n,1);
            for k=1:n
                thisyv{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
            end
            ylabelv = glist;

            fw = gui.myWaitbar(FigureHandle);
            % gui.sc_uitabgrpfig_expplot(thisyv, ylabelv, sce.s, FigureHandle);
            gui.sc_uitabgrpfig_feaplot(thisyv, ylabelv, sce.s, FigureHandle, 2);   
            gui.myWaitbar(FigureHandle, fw);
            
        case 'Cell State'
            [thisyv, ylabelv] = gui.i_selectnstates(sce, true, [], FigureHandle);

            a = false(length(thisyv), 1);
            for k = 1:length(thisyv)
                a(k) = isnumeric(thisyv{k});
            end
            if any(a)
                if ~all(a)
                    thisyv = thisyv(a);
                    ylabelv = ylabelv(a);
                    gui.myHelpdlg(FigureHandle, 'Only continuous variables of cell state will be shown.');
                end                
            else
                gui.myHelpdlg(FigureHandle, 'No valid cell state variables.');
                return;
            end
            gui.sc_uitabgrpfig_feaplot(thisyv, ylabelv, sce.s, FigureHandle, 2);    

        otherwise
            return;
    end

end
