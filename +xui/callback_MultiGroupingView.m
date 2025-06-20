function callback_MultiGroupingView(app, ~)
    answer = gui.myQuestdlg(app.UIFigure, 'Select type of multi-view:','', ...
        {'Multigrouping','Multiembedding'},'Multigrouping');

    switch answer
        case 'Multigrouping'
            if gui.i_isuifig(app.UIFigure), focus(app.UIFigure); end
            [thiscv, clabelv] = gui.i_selectnstates(app.sce, false, [4, 5], app.UIFigure);
            if isempty(thiscv) || isempty(clabelv), return; end

            hx=gui.myFigure;
            hFig = hx.FigHandle;
            hFig.Position(3) = hFig.Position(3) * 1.8;
            axesv = cell(length(thiscv),1);
            cmapv = cell(length(thiscv),1);
            hv = cell(length(thiscv),1);

            for k = 1:length(thiscv)
                axesv{k} = nexttile;
                hv{k} = gui.i_gscatter3(app.sce.s, thiscv{k}, 1, 1);
                title(strrep(clabelv{k},'_','\_'));
                cmapv{k} = colormap;
            end

            hx.addCustomButton('off', @in_showclustlabel, "label.jpg", "Show cluster labels");
            hx.show(app.UIFigure);

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
            listitems = fieldnames(app.sce.struct_cell_embeddings);
            n = length(listitems);
            valididx = false(n,1);
            for k=1:n
                s = app.sce.struct_cell_embeddings.(listitems{k});
                if ~isempty(s) && size(s,2)>1 && size(s,1)==app.sce.NumCells
                    valididx(k)=true;
                end
            end
            listitems = listitems(valididx);
            if isempty(listitems)
                gui.myWarndlg(app.UIFigure, 'No embeding is available.');
                return;
            end
            n = length(listitems);

            if gui.i_isuifig(app.UIFigure)
                [indx2, tf2] = gui.myListdlg(app.UIFigure, listitems, ...
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
                gui.sc_multiembeddingview(app.sce, listitems(indx2), app.UIFigure);
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

    function in_showclustlabel(~, ~)
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
                    siv = app.sce.s(idx, :);
                    si = mean(siv, 1);
                    [kk] = dsearchn(siv, si);
                    datatip(hv{kx}, 'DataIndex', idx(kk));
                end        
            end
        end
    end

end
