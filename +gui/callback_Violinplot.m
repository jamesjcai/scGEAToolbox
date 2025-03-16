function callback_Violinplot(src, ~)

    if isa(src,"SingleCellExperiment")
        sce = src;
        FigureHandle = [];
    else
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end
    
    [thisc, ~] = gui.i_select1class(sce,[],[],[],FigureHandle);
    if isempty(thisc), return; end


    [~, cLorder] = grp2idx(thisc);
    [newidx] = gui.i_selmultidlg(cLorder, cLorder, FigureHandle);
    if isempty(newidx), return; end
    picked = ismember(thisc, cLorder(newidx));
    % cLorderx = cLorder(ismember(cLorder,cLorder(newidx)));
    if ~all(picked), thisc = thisc(picked); end
    
    answer = gui.myQuestdlg(FigureHandle, "Violinplot for gene expression or cell state variables?","", ...
        {'Gene Expression', 'Cell State'},'Gene Expression');
    switch answer
        case 'Gene Expression'
            % [c, cL] = grp2idx(thisc);
            % [c, cL, noanswer] = gui.i_reordergroups(thisc);
            % if noanswer, return; end
            [glist] = gui.i_selectngenes(sce, [], FigureHandle);
            if isempty(glist)
                gui.myHelpdlg(FigureHandle, 'No gene selected.');
                return;
            end

            [Xt] = gui.i_transformx(sce.X, true, 8, FigureHandle);
            if isempty(Xt), return; end

            if isscalar(glist)
                answer='No';
            else
                answer = gui.myQuestdlg(FigureHandle, 'Plot all in the same figure?','');
            end

            if strcmp(answer, 'Yes')
                [c, cL, noanswer] = gui.i_reordergroups(thisc, [], FigureHandle);
                if noanswer, return; end
                % [c, cL] = grp2idx(thisc);
                fw = gui.myWaitbar(FigureHandle);
                gui.i_violinmatrix(full(Xt), sce.g, c, cL, glist, ...
                        false, '', FigureHandle);
                gui.myWaitbar(FigureHandle, fw);

                return;
            end

            n = length(glist);
            thisyv = cell(n,1);
            for k=1:n
                thisyv{k} = full(Xt(upper(sce.g) == upper(glist(k)), :));
                if ~all(picked)
                    thisyv{k} = thisyv{k}(picked);
                end
            end
            ylabelv = glist;

        case 'Cell State'
            [thisyv, ylabelv] = gui.i_selectnstates(sce, true);

            a = false(length(thisyv), 1);
            for k = 1:length(thisyv)
                thisyv{k} = thisyv{k}(picked);
                a(k) = isnumeric(thisyv{k});
            end
            if any(a)
                if ~all(a)
                    thisyv = thisyv(a);
                    ylabelv = ylabelv(a);
                    gui.myHelpdlg(FigureHandle, 'Only continuous variables of cell state will be shown.');
                end                
            else
                gui.myHelpdlg(FigureHandle, 'No valid cell state variables. Violinplot cannot be shown.');
            end
        otherwise
            return;
    end
    gui.sc_uitabgrpfig_vioplot(thisyv, ylabelv, thisc, FigureHandle);
end