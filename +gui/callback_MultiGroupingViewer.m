function callback_MultiGroupingViewer(src, ~)
    FigureHandle = src.Parent.Parent;
    sce = guidata(FigureHandle);
    
    answer = questdlg('Select type of multi-view:','', ...
        'Multigrouping','Multiembedding','Multigrouping');

    switch answer
        % case 'Two-group'
        % 
        %     if matlab.ui.internal.isUIFigure(FigureHandle), focus(FigureHandle); end
        %     [thisc1, clabel1, thisc2, clabel2] = gui.i_select2state_new(sce);
        %     if isempty(thisc1) || isempty(thisc2), return; end
        % 
        %     if matlab.ui.internal.isUIFigure(FigureHandle), focus(FigureHandle); end    
        %     fw=gui.gui_waitbar;    
        %     [c, cL] = grp2idx(thisc1);
        %     cx1.c = c;
        %     cx1.cL = strrep(cL, '_', '\_');
        %     [c, cL] = grp2idx(thisc2);
        %     cx2.c = c;
        %     cx2.cL = strrep(cL, '_', '\_');
        %     gui.sc_multigroupings(sce, cx1, cx2, clabel1, clabel2, FigureHandle);
        %     gui.gui_waitbar(fw);
        % 
        %     if matlab.ui.internal.isUIFigure(FigureHandle), focus(FigureHandle); end

        case 'Multigrouping'
            if matlab.ui.internal.isUIFigure(FigureHandle), focus(FigureHandle); end
            [thiscv, clabelv] = gui.i_selectnstates(sce);
            if isempty(thiscv) || isempty(clabelv), return; end

            hFig = figure('Visible','off');
            hFig.Position(3) = hFig.Position(3) * 1.8;
            axesv = cell(length(thiscv),1);
            for k=1:length(thiscv)
                axesv{k} = nexttile;
                gui.i_gscatter3(sce.s, thiscv{k}, 1, 1);
                title(strrep(clabelv{k},'_','\_'));
            end

            gui.i_movegui2parent(hFig, FigureHandle);

            drawnow;
            hFig.Visible=true;
            dt = datacursormode(hFig);
            dt.UpdateFcn = {@in_myupdatefcnx12};
            %evalin('base', 'h = findobj(gcf,''type'',''axes'');');
            %evalin('base', 'hlink = linkprop(h, {''CameraPosition'',''CameraUpVector''});');
            evalin('base', 'linkprop(findobj(gcf,''type'',''axes''), {''CameraPosition'',''CameraUpVector''});');
            %h = findobj(hFig,'type','axes');
            %linkprop(h, {'CameraPosition','CameraUpVector'});
            rotate3d(hFig,'on');
            hBr = brush(hFig);
            hBr.ActionPostCallback = {@onBrushAction, axesv};

            tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
            uipushtool(tb, 'Separator', 'off');
            gui.gui_3dcamera(tb, 'AllCells');
            
        case 'Multiembedding'
            listitems = fieldnames(sce.struct_cell_embeddings);
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
                'InitialValue', 1:n, ...
                'ListSize', [220, 300]);
            if tf2 == 1
                gui.sc_multiembeddingview(sce, listitems(indx2), FigureHandle);
            end
        otherwise

    end

        function onBrushAction(~, event, axv)
            for kx=1:length(axv)
                if isequal(event.Axes, axv{kx})
                    idx = kx;
                    continue;
                end
            end
            d = axv{idx}.Children.BrushData;
            for kx=1:length(axv)
                if kx ~= idx
                    axv{kx}.Children.BrushData = d;
                end
            end
        end

    function [txt] = in_myupdatefcnx12(Targxet, event_obj)
        % pos = event_obj.Position;
       

        for kx=1:length(axesv)
            if isequal(Targxet.Parent, axesv{kx})
                idx = event_obj.DataIndex;
                c1 = thiscv{kx};
                %[~,cL1]=grp2idx(c1);                   
                txt = c1(idx);
                if isstring(txt) || ischar(txt)
                    txt = strrep(txt,'_','\_');
                end
                continue;
            end
        end
    end    
end
