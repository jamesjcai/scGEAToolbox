function [requirerefresh] = callback_RenameCellType(src)
requirerefresh = false;

[FigureHandle, sce] = gui.gui_getfigsce(src);

if isempty(sce.c_cell_type_tx)
    gui.myErrordlg(FigureHandle, 'sce.c_cell_type_tx undefined');
    return;
end

removesubscript = false;
if any(contains(sce.c_cell_type_tx, "_{"))
    answer = gui.myQuestdlg(FigureHandle, ...
        'Remove subscript of cell type names?', ...
        '', {'Yes', 'No', 'Cancel'}, 'No');
    switch answer
        case 'Yes'
            removesubscript = true;
        case 'No'
            removesubscript = false;
        otherwise
            return;
    end
end
if removesubscript

    gui.myHelpdlg(FigureHandle, 'Subscript will be removed.');
    if ~isstring(sce.c_cell_type_tx)
        sce.c_cell_type_tx = string(sce.c_cell_type_tx);
    end
    for k = 1:length(sce.c_cell_type_tx)
        a = strfind(sce.c_cell_type_tx(k), '_{');
        if ~isempty(a)
            sce.c_cell_type_tx(k) = extractBefore(sce.c_cell_type_tx(k), a);
        end
    end
    requirerefresh = true;
else
    [ci, cLi] = pkg.i_grp2idxsorted(sce.c_cell_type_tx);

    % %---------
    % [cLisorted,idx]=natsort(cLi);
    % cisorted=ci;
    % for k=1:length(idx), cisorted(ci==idx(k))=k; end
    % ci=cisorted;
    % cLi=cLisorted;
    % %----------
    if gui.i_isuifig(FigureHandle)
        [indxx, tfx] = gui.myListdlg(FigureHandle, cLi, ...
            'Select cell type');
    else
        [indxx, tfx] = listdlg('PromptString', ...
            {'Select cell type'}, ...
            'SelectionMode', 'single', ...
            'ListString', string(cLi), 'ListSize', [220, 300]);
    end
    if tfx == 1
        i = ismember(ci, indxx);
        
if gui.i_isuifig(FigureHandle)
    newctype = gui.myInputdlg({'New cell type'}, 'Rename', ...
        cLi(ci(i)), FigureHandle);
else
    newctype = inputdlg('New cell type', 'Rename', [1, 50], cLi(ci(i)));
end


        if ~isempty(newctype)
            cLi(ci(i)) = newctype;
            sce.c_cell_type_tx = string(cLi(ci));
            requirerefresh = true;
            %[c, cL] = grp2idx(sce.c_cell_type_tx);
            %i_labelclusters(false);
        end
    end
end
guidata(FigureHandle, sce);
end
