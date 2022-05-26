function [requirerefresh]=callback_RenameCellType(src)
    requirerefresh=false;
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);

    if isempty(sce.c_cell_type_tx)
        errordlg('sce.c_cell_type_tx undefined');
        return;
    end
    %answer = questdlg('Rename a cell type?');
    %if ~strcmp(answer, 'Yes'), return; end
    [ci, cLi] = grp2idx(sce.c_cell_type_tx);
    [indxx, tfx] = listdlg('PromptString',...
        {'Select cell type'},...
        'SelectionMode', 'single',...
        'ListString', string(cLi));
    if tfx == 1
        i = ismember(ci, indxx);
        newctype = inputdlg('New cell type', 'Rename', [1 50], cLi(ci(i)));
        if ~isempty(newctype)
            cLi(ci(i)) = newctype;
            sce.c_cell_type_tx = string(cLi(ci));
            requirerefresh=true;
            %[c, cL] = grp2idx(sce.c_cell_type_tx);
            %i_labelclusters(false);
        end
    end
    guidata(FigureHandle, sce);
end
