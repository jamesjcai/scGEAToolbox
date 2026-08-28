function callback_ShowGeneExprCompr(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

[axx, bxx] = view(findall(FigureHandle,'type','axes'));

[glist] = gui.i_selectngenes(sce, [], FigureHandle);
if isempty(glist), return; end

allowunique = false;
% Several grouping variables may be picked; they cross into one composite
% label per cell ("Macrophages | IL"). Downstream treats thisc as a
% per-cell label vector, so the composite needs no special handling.
[thisc] = gui.i_selectnclass(sce, allowunique,[],[],FigureHandle);
if isempty(thisc), return; end
if isscalar(unique(thisc))
        answer = gui.myQuestdlg(FigureHandle, "All cells are in the same group. No comparison will be made. Continue?", ...
            "", {'Yes', 'No', 'Cancel'}, 'No');
        switch answer
            case 'Yes'
            otherwise
                return;
        end
    else    % length(unique(thisc)) ~= 1
        [ci, cLi] = findgroups(string(thisc));
        listitems = natsort(cLi);
        n = length(listitems);

    if gui.i_isuifig(FigureHandle)
        [indxx, tfx] = gui.myListdlg(FigureHandle, ...
            listitems, 'Select two groups:', ...
            listitems);
    else
        [indxx, tfx] = listdlg('PromptString', ...
            {'Select two groups:'}, ...
            'SelectionMode', 'multiple', ...
            'ListString', listitems, ...
            'InitialValue', 1:n, 'ListSize', [220, 300]);
    end

        if tfx == 1
            [y1, idx1] = ismember(listitems(indxx), cLi);
            assert(all(y1));
            idx2 = ismember(ci, idx1);
            sce1 = copy(sce).selectcells(idx2);  % OK
            thisc = thisc(idx2);
        else
            return;
        end
    end


fw=gui.myWaitbar(FigureHandle);
gui.sc_uitabgrpfig_expcomp(sce1, glist, FigureHandle, [axx, bxx], thisc);
gui.myWaitbar(FigureHandle, fw);

end
