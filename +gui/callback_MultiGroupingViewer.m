function callback_MultiGroupingViewer(src, ~)
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    
    answer = questdlg('Select the type of multi-view:','', ...
        'Multigrouping','Multiembedding','Cancel','Multigrouping');

    switch answer
        case 'Multigrouping'

            if matlab.ui.internal.isUIFigure(FigureHandle), focus(FigureHandle); end
            [thisc1, clable1, thisc2, clable2] = gui.i_select2state_new(sce);
            if isempty(thisc1) || isempty(thisc2), return; end
            
            if matlab.ui.internal.isUIFigure(FigureHandle), focus(FigureHandle); end    
            fw=gui.gui_waitbar;    
            [c, cL] = grp2idx(thisc1);
            cx1.c = c;
            cx1.cL = strrep(cL, '_', '\_');
            [c, cL] = grp2idx(thisc2);
            cx2.c = c;
            cx2.cL = strrep(cL, '_', '\_');
            gui.sc_multigroupings(sce, cx1, cx2, clable1, clable2, FigureHandle);
            gui.gui_waitbar(fw);
            
            if matlab.ui.internal.isUIFigure(FigureHandle), focus(FigureHandle); end

        case 'Multiembedding'
            listitems = fields(sce.struct_cell_embeddings);
            n = length(listitems);
            valididx = false(n,1);
            for k=1:n
                s = sce.struct_cell_embeddings.(listitems{k});
                if ~isempty(s) && size(s,2)>1 && size(s,1)==sce.NumCells
                    valididx(k)=true;
                end
            end
            listitems = listitems(valididx);
            if isempty(listitems)
                warndlg('No embeding is available.','');
                return;
            end
            n = length(listitems);

            [indx2, tf2] = listdlg('PromptString', ...
                {'Select embeddings:'}, ...
                'SelectionMode', 'multiple', ...
                'ListString', listitems, ...
                'InitialValue', 1:n);
            if tf2 == 1
                hFig = figure('Visible', false);                
                hFig.Position(3) = hFig.Position(3) * 1.8;
                [~, newpos] = gui.i_getchildpos(hFig, FigureHandle);
                movegui(hFig, newpos);
                figure(hFig);
                for k=1:length(indx2)
                    s = sce.struct_cell_embeddings.(listitems{k});
                    if size(s,2)>1 && size(s,1)==sce.NumCells
                        nexttile
                        gui.i_gscatter3(s, sce.c, 1, 1);
                        title(listitems{k});
                    end
                end
                hFig.Visible=true;
            end            
        otherwise
    end
end
