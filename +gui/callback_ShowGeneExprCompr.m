function callback_ShowGeneExprCompr(src, ~)


    [FigureHandle, sce, isui] = gui.gui_getfigsce(src);
    [axx, bxx] = view(findall(FigureHandle,'type','axes'));
    
    [glist] = gui.i_selectngenes(sce, [], FigureHandle);
    if isempty(glist), return; end

        allowunique = false;
        [thisc] = gui.i_select1class(sce, allowunique);
        if isempty(thisc), return; end
        if isscalar(unique(thisc))
            answer = questdlg("All cells are in the same group. No comparison will be made. Continue?", ...
                "", 'Yes', 'No', 'Cancel', 'No');
            switch answer
                case 'Yes'
                otherwise
                    return;
            end
        else    % length(unique(thisc)) ~= 1
            [ci, cLi] = grp2idx(thisc);
            listitems = natsort(string(cLi));
            n = length(listitems);
            [indxx, tfx] = listdlg('PromptString', {'Select two groups:'}, ...
                'SelectionMode', 'multiple', ...
                'ListString', listitems, ...
                'InitialValue', 1:n, 'ListSize', [220, 300]);
            if tfx == 1
                [y1, idx1] = ismember(listitems(indxx), cLi);
                assert(all(y1));
                idx2 = ismember(ci, idx1);
                sce = sce.selectcells(idx2);
                thisc = thisc(idx2);
            else
                return;
            end
        end


    fw=gui.gui_waitbar;
    gui.sc_uitabgrpfig_expcomp(sce, glist, FigureHandle, [axx, bxx], thisc);
    gui.gui_waitbar(fw);

end
