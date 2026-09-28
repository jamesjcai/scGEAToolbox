function sc_multiembeddingview(sce, embeddingtags, parentfig)

if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
if isempty(embeddingtags)
    embeddingtags = fieldnames(sce.struct_cell_embeddings);
end
hx=gui.myFigure(parentfig);
hFig=hx.FigHandle;
hFig.Position(3) = hFig.Position(3) * 1.8;
axesv = cell(length(embeddingtags),1);
for k = 1:length(embeddingtags)
        s = sce.struct_cell_embeddings.(embeddingtags{k});
        if size(s,2)>1 && size(s,1)==sce.NumCells
            axesv{k} = nexttile;
            gui.i_gscatter3(s, sce.c, 1, 1, axesv{k});
            title(axesv{k}, embeddingtags{k});
        end
    end


hBr = brush(hFig);
hBr.ActionPostCallback = {@onBrushAction, axesv};

hx.addCustomButton('off',  @in_showgeneexp, 'google-docs.jpg', 'Select a gene to show expression...');
hx.addCustomButton('off',  @in_showcellstate, 'bookmark-book.jpg', 'Show cell state...');
hx.show(parentfig);


function in_showcellstate(~, ~)
        [thisc, clabel] = gui.i_select1state(sce, false, false, true, false, hFig);
        if isempty(thisc), return; end
        [c, cL] = findgroups(string(thisc));
        stxtyes = cL(c);
        if isstring(stxtyes) || iscellstr(stxtyes)
            stxtyes = strrep(stxtyes, "_", "\_");
            stxtyes = strtrim(stxtyes);
        end
        row = dataTipTextRow('', stxtyes);
        a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
        for kx = 1:length(axesv)
           s = sce.struct_cell_embeddings.(embeddingtags{kx});
           if ~isempty(s) && ~isempty(axesv{kx})
               h = gui.i_gscatter3(s, c, 1, 1, axesv{kx});
               title(axesv{kx}, string(embeddingtags{kx})+" - "+string(clabel));
               h.DataTipTemplate.DataTipRows = row;
               gui.i_setautumncolor(c, a, true, false, axesv{kx}, hFig);
           end
        end
    end

function in_showgeneexp(~, ~)
        [gsorted] = gui.i_sortgenenames(sce, hFig);
        if isempty(gsorted), return; end
        figure(hFig);
       if gui.i_isuifig(hFig)
            [indx, tf] = gui.myListdlg(hFig, gsorted, 'Select a gene:', [], false);
        else
            [indx, tf] = listdlg('PromptString', 'Select a gene:', ...
                'SelectionMode', 'single', 'ListString', ...
                gsorted, 'ListSize', [220, 300]);
        end

        if tf == 1
            c = full(sce.X(sce.g == gsorted(indx), :));
            % Its own: A was only ever set inside in_showcellstate, where it
            % is local, so this button always errored on the first panel.
            a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

            for kx = 1:length(axesv)
               if isempty(axesv{kx}), continue; end   % embedding not drawn
               s = sce.struct_cell_embeddings.(embeddingtags{kx});
               gui.i_gscatter3(s, c, 1, 1, axesv{kx});
               title(axesv{kx}, string(embeddingtags{kx})+" - "+string(gsorted(indx)));
               gui.i_setautumncolor(c, a, true, any(c==0), axesv{kx}, hFig);
            end
        end
    end

function onBrushAction(~, event, axv)
        % Copy the brushed cells to every other panel. By the panel's
        % scatter, not .Children: that is several objects once a panel
        % holds anything else, and IDX was undefined for an axes not in AXV.
        src = [];
        for kx = 1:length(axv)
            if ~isempty(axv{kx}) && isequal(event.Axes, axv{kx})
                src = findobj(axv{kx}, 'Type', 'scatter');
                break;
            end
        end
        if isempty(src), return; end
        d = src(1).BrushData;
        for kx = 1:length(axv)
            if isempty(axv{kx}) || isequal(event.Axes, axv{kx}), continue; end
            h = findobj(axv{kx}, 'Type', 'scatter');
            if ~isempty(h) && numel(h(1).BrushData) == numel(d)
                h(1).BrushData = d;
            end
        end
    end

end
