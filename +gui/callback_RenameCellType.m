function [requirerefresh] = callback_RenameCellType(src)
requirerefresh = false;
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

if isempty(sce.c_cell_type_tx)
    errordlg('sce.c_cell_type_tx undefined');
    return;
end

removesubscript = false;
if any(contains(sce.c_cell_type_tx, "_{"))
    answer = questdlg('Remove subscript of cell type names?', ...
        '', 'Yes', 'No', 'Cancel', 'No');
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

    waitfor(helpdlg('Subscript will be removed.'));
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

    [indxx, tfx] = listdlg('PromptString', ...
        {'Select cell type'}, ...
        'SelectionMode', 'single', ...
        'ListString', string(cLi), 'ListSize', [220, 300]);
    if tfx == 1
        i = ismember(ci, indxx);
        newctype = inputdlg('New cell type', 'Rename', ...
            [1, 60], cLi(ci(i)));
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
