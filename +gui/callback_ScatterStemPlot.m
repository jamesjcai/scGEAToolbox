function callback_ScatterStemPlot(src, ~)

    if isa(src,"SingleCellExperiment")
        sce = src;
        FigureHandle = [];
    else
        [FigureHandle, sce, isui] = gui.gui_getfigsce(src);
    end
    
    %[thisc, ~] = gui.i_select1class(sce);
    %if isempty(thisc), return; end
    %[~, cLorder] = grp2idx(thisc);
    %[newidx] = gui.i_selmultidlg(cLorder, cLorder, FigureHandle);
    %if isempty(newidx), return; end
    %picked = ismember(thisc, cLorder(newidx));
    %if ~all(picked), thisc = thisc(picked); end
    
    answer = questdlg("Scatter-stem plot for gene expression or cell state variables?","", ...
        'Gene Expression', 'Cell State', 'Gene Expression');
    switch answer
        case 'Gene Expression'
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist)
                helpdlg('No gene selected.', '');
                return;
            end

            [Xt] = gui.i_transformx(sce.X, [], [], FigureHandle);
            if isempty(Xt), return; end

            % answer = questdlg('Plot all in the same figure?','');
            % if strcmp(answer, 'Yes')                
            %     fw = gui.gui_waitbar;                
            %     gui.i_violinmatrix(full(Xt), sce.g, c, cL, glist, ...
            %             false, '', FigureHandle);
            % 
            %     gui.gui_waitbar(fw);
            % 
            %     return;
            % end

            n = length(glist);
            thisyv = cell(n,1);
            for k=1:n
                thisyv{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
            end
            ylabelv = glist;

        case 'Cell State'
            [thisyv, ylabelv] = gui.i_selectnstates(sce, true);

            a = false(length(thisyv), 1);
            for k = 1:length(thisyv)
                a(k) = isnumeric(thisyv{k});
            end
            if any(a)
                if ~all(a)
                    thisyv = thisyv(a);
                    ylabelv = ylabelv(a);
                    waitfor(helpdlg('Only continuous variables of cell state will be shown.',''));
                end                
            else
                waitfor(helpdlg('No valid cell state variables. Violinplot cannot be shown.',''));
            end
        otherwise
            return;
    end
    gui.sc_uitabgrpfig_feaplot(thisyv, ylabelv, sce.s, FigureHandle, 2);    
end
