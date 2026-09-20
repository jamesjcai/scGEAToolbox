function callback_MultiGroupingView(src, ~)


[FigureHandle, sce] = gui.gui_getfigsce(src);
answer = gui.myQuestdlg(FigureHandle, 'Select type of multi-view:','', ...
{'Multigrouping','Multiembedding'},'Multigrouping');

% figure(FigureHandle);
gui.i_bringtofront(FigureHandle);
switch answer
    case 'Multigrouping'
        preftagname ='selectednstates';
        defaultindx = getpref('scgeatoolbox', preftagname, [4, 5]);
        [thiscv, clabelv] = gui.i_selectnstates(sce, false, ...
            defaultindx, FigureHandle);
        if isempty(thiscv) || isempty(clabelv), return; end

        % The panels themselves live in GUI.I_MULTIGROUPVIEW, so that a
        % caller that already knows which groupings it wants - the cell type
        % annotation comparison - draws the same linked figure without this
        % picker in front of it.
        gui.i_multigroupview(sce, thiscv, string(clabelv), FigureHandle, ...
            'Multi-Grouping View');

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
            figure(FigureHandle);
            gui.sc_multiembeddingview(sce, listitems(indx2), FigureHandle);
        end
    otherwise
        return;
end


end
