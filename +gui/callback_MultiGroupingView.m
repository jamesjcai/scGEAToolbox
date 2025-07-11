function callback_MultiGroupingView(src, ~)
    

        [FigureHandle, sce] = gui.gui_getfigsce(src);
    answer = gui.myQuestdlg(FigureHandle, 'Select type of multi-view:','', ...
        {'Multigrouping','Multiembedding'},'Multigrouping');

    switch answer
        % case 'Two-group'
        % 
        %     [thisc1, clabel1, thisc2, clabel2] = gui.i_select2states_new(sce);
        %     if isempty(thisc1) || isempty(thisc2), return; end
        %         
        %     fw=gui.myWaitbar(FigureHandle);    
        %     [c, cL] = grp2idx(thisc1);
        %     cx1.c = c;
        %     cx1.cL = strrep(cL, '_', '\_');
        %     [c, cL] = grp2idx(thisc2);
        %     cx2.c = c;
        %     cx2.cL = strrep(cL, '_', '\_');
        %     gui.sc_multigroupings(sce, cx1, cx2, clabel1, clabel2, FigureHandle);
        %     gui.myWaitbar(FigureHandle, fw);
        % 
        

        case 'Multigrouping'
            if gui.i_isuifig(FigureHandle), focus(FigureHandle); end
            [thiscv, clabelv] = gui.i_selectnstates(sce, false, [4, 5], FigureHandle);
            if isempty(thiscv) || isempty(clabelv), return; end

            hx=gui.myFigure(FigureHandle);
            hFig = hx.FigHandle;
            hFig.Position(3) = hFig.Position(3) * 1.8;
            axesv = cell(length(thiscv),1);
            cmapv = cell(length(thiscv),1);
            hv = cell(length(thiscv),1);

            for k = 1:length(thiscv)
                axesv{k} = nexttile;
                hv{k} = gui.i_gscatter3(sce.s, thiscv{k}, 1, 1);
                title(strrep(clabelv{k},'_','\_'));
                cmapv{k} = colormap;
            end

            hx.addCustomButton('off', @in_callback_showclustlabel, "label.jpg", "Show cluster labels");
            hx.show(FigureHandle);

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

            for k = 1:length(thiscv)
               colormap(axesv{k}, cmapv{k});
            end            
            
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
                gui.myWarndlg(FigureHandle, 'No embeding is available.');
                return;
            end
            n = length(listitems);

            if gui.i_isuifig(FigureHandle)
                [indx2, tf2] = gui.myListdlg(FigureHandle, listitems, ...
                    'Select embeddings:', ...
                    listitems);
            else            
                [indx2, tf2] = listdlg('PromptString', ...
                    {'Select embeddings:'}, ...
                    'SelectionMode', 'multiple', ...
                    'ListString', listitems, ...
                    'InitialValue', 1:n, ...
                    'ListSize', [220, 300]);
            end
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

    function in_callback_showclustlabel(~, ~)
        hastip = false;
        for kx = 1:length(thiscv)
            dtp1 = findobj(hv{kx}, 'Type', 'datatip');
            if ~isempty(dtp1)
                delete(dtp1);
                hastip = true;                
            end
        end
        if hastip, return; end

        for kx = 1:length(thiscv)
            [c1, cL1] = grp2idx(thiscv{kx});
            cL1 = strrep(cL1,'_','\_');
            if max(c1) < 50
                hv{kx}.DataTipTemplate.DataTipRows = dataTipTextRow('', cL1(c1));
                for i = 1:max(c1)
                    idx = find(c1 == i);
                    siv = sce.s(idx, :);
                    si = mean(siv, 1);
                    [kk] = dsearchn(siv, si);
                    datatip(hv{kx}, 'DataIndex', idx(kk));
                end        
            end
        end
    end

end
