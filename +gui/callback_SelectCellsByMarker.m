function callback_SelectCellsByMarker(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

answer1 = questdlg('Use single or mulitple markers?', ...
    'Single/Multiple Markers', 'Single', 'Multiple', 'Cancel', 'Single');
switch answer1
    case 'Single'
        do_single;
    case 'Multiple'
        do_multiple;
    otherwise
        return;
end

    function do_multiple
        [glist] = gui.i_selectngenes(sce);
        if ~isempty(glist)
            [y, i] = ismember(upper(glist), upper(sce.g));
            if ~all(y), error('Unspecific running error.'); end
            ix = sum(sce.X(i, :) > 0, 1) == length(i);
            if ~any(ix)
                helpdlg('No cells expressing all selected markers.', '');
                return;
            end

            answer = questdlg('Extract markers+ or markers- cells?', ...
                'Positive or Negative', ...
                'Markers+', 'Markers-', 'Cancel', 'Markers+');
            switch answer
                case 'Markers+'
                    idx = ix;
                case 'Markers-'
                    idx = ~ix;
                case 'Cancel'
                    return;
                otherwise
                    return;
            end
            [ax, bx] = view();
            fw = gui.gui_waitbar;
            scex = selectcells(sce, idx);
            sc_scatter_sce(scex);
            view(ax, bx);
            gui.gui_waitbar(fw);

        end
    end


    function do_single
        [gsorted] = gui.i_sortgenenames(sce);
        if isempty(gsorted), return; end

        [indx, tf] = listdlg('PromptString', {'Select a gene', '', ''}, ...
            'SelectionMode', 'single', 'ListString', gsorted);
        if tf == 1
            [ax, bx] = view();
            tg = gsorted(indx);
            c = sce.X(sce.g == tg, :);
            answer = questdlg(sprintf('Extract %s+ or %s- cells?', tg, tg), 'Positive or Negative', ...
                sprintf('%s+', tg), sprintf('%s-', tg), 'Split', sprintf('%s+', tg));
            if strcmp(answer, sprintf('%s+', tg))
                idx = c > 0;
            elseif strcmp(answer, sprintf('%s-', tg))
                idx = c == 0;
            elseif strcmp(answer, 'Split')
                idx1 = c > 0;
                idx2 = c == 0;

                scex = sce;
                scex.list_cell_attributes = [sce.list_cell_attributes, {'old_batch_id', scex.c_batch_id}];
                scex.c_batch_id(idx1) = sprintf('%s+', tg);
                scex.c_batch_id(idx2) = sprintf('%s-', tg);
                scex.c = scex.c_batch_id;
                sc_scatter_sce(scex);

                %{
                scex = selectcells(sce, idx1);
                fx = sc_scatter_sce(scex);
                fx.Position(3:4) = 0.8 * fx.Position(3:4);
                movegui(fx, 'center');
                fx.Position(1) = fx.Position(1) - 250;
                fx = fx.CurrentAxes;
                fx.Subtitle.String = sprintf('%s\n%s', fx.Subtitle.String, sprintf('%s+', tg));
                view(fx, ax, bx);

                answer = questdlg(sprintf('%s Cells extracted. Continue?', sprintf('%s+', tg)), '');
                if ~strcmp(answer, 'Yes'), return; end
                scey = selectcells(sce, idx2);
                fy = sc_scatter_sce(scey);
                fy.Position(3:4) = 0.8 * fy.Position(3:4);
                movegui(fy, 'center');
                fy.Position(1) = fy.Position(1) + 250;                
                fy = fy.CurrentAxes;
                fy.Subtitle.String = sprintf('%s\n%s', fy.Subtitle.String, sprintf('%s-', tg));
                view(fy, ax, bx);
                %}

                return;                
            elseif strcmp(answer, 'Cancel')
                return;
            else
                return;
            end
            fw = gui.gui_waitbar;
            scex = selectcells(sce, idx);
            sc_scatter_sce(scex);
            view(ax, bx);
            gui.gui_waitbar(fw);
        end
    end
end
