function callback_SelectCellsByMarker(src, ~)
% The extracted cells replace the dataset in this window; the menu handler
% takes the Undo snapshot, so Ctrl+Z brings the rest back. They used to open
% in a second scgeatoolApp window, over this one.

[FigureHandle, sce] = gui.gui_getfigsce(src);

answer1 = gui.myQuestdlg(FigureHandle, 'Use single or mulitple markers?', ...
'Single/Multiple Markers', {'Single', 'Multiple', 'Cancel'}, 'Single');
switch answer1
    case 'Single'
        do_single;
    case 'Multiple'
        do_multiple;
    otherwise
        return;
end

function do_multiple
        [glist] = gui.i_selectngenes(sce, [] , FigureHandle);
        if ~isempty(glist)
            [y, i] = ismember(upper(glist), upper(sce.g));
            if ~all(y), error('Unspecific running error.'); end
            ix = sum(sce.X(i, :) > 0, 1) == length(i);
            if ~any(ix)
                gui.myHelpdlg(FigureHandle, 'No cells expressing all selected markers.', '');
                return;
            end

            answer = gui.myQuestdlg(FigureHandle, 'Extract markers+ or markers- cells?', ...
                'Positive or Negative', ...
                {'Markers+', 'Markers-', 'Cancel'}, 'Markers+');
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
            fw = gui.myWaitbar(FigureHandle);
            scex = copy(sce).selectcells(idx);  % OK
            gui.myWaitbar(FigureHandle, fw);
            gui.i_replacesce(src, scex, sce.NumCells);

        end
    end


function do_single
        [gsorted] = gui.i_sortgenenames(sce, FigureHandle);
        if isempty(gsorted), return; end

       if gui.i_isuifig(FigureHandle)
           [indx, tf] = gui.myListdlg(FigureHandle, gsorted, ...
                'Select a gene', [], false);
        else
            [indx, tf] = listdlg('PromptString', {'Select a gene', '', ''}, ...
                'SelectionMode', 'single', ...
                'ListString', gsorted, 'ListSize', [220, 300]);
       end
        if tf == 1
            tg = gsorted(indx);
            c = sce.X(sce.g == tg, :);
            answer = gui.myQuestdlg(FigureHandle, sprintf('Extract %s+ or %s- cells?', tg, tg), ...
                'Positive or Negative', ...
                {sprintf('%s+', tg), sprintf('%s-', tg), 'Split'}, sprintf('%s+', tg));
            if strcmp(answer, sprintf('%s+', tg))
                idx = c > 0;
            elseif strcmp(answer, sprintf('%s-', tg))
                idx = c == 0;
            elseif strcmp(answer, 'Split')
                idx1 = c > 0;
                idx2 = c == 0;

                scex = copy(sce);
                scex.setCellAttribute('old_batch_id', scex.c_batch_id);
                scex.c_batch_id(idx1) = sprintf('%s+', tg);
                scex.c_batch_id(idx2) = sprintf('%s-', tg);
                scex.c = scex.c_batch_id;
                % Every cell is kept, only relabelled, so there is nothing
                % to report the way a subset reports its cell count.
                src.sce = scex;
                [src.c, src.cL] = findgroups(string(scex.c));
                src.in_RefreshAll(true, false);
                return;
            elseif strcmp(answer, 'Cancel')
                return;
            else
                return;
            end
            fw = gui.myWaitbar(FigureHandle);
            scex = copy(sce).selectcells(idx); % OK
            gui.myWaitbar(FigureHandle, fw);
            gui.i_replacesce(src, scex, sce.NumCells);
        end
    end
end
