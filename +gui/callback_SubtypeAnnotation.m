function [requirerefresh] = callback_SubtypeAnnotation(src, ~)

    requirerefresh = false;
    % answer = gui.myQuestdlg(FigureHandle, 'This function assigns subtype name to selected cell type, continue?', '');
    % if ~strcmp(answer, 'Yes'), return; end

    if isa(src, "SingleCellExperiment")
        sce = src;
    else
        [FigureHandle, sce] = gui.gui_getfigsce(src);
    end

    pw1 = fileparts(mfilename('fullpath'));

    pth2 = fullfile(pw1, '..', 'resources', 'PanglaoDB', 'cellsubtypes.xlsx');
    T = readtable(pth2);
    ctypes = unique(string(sce.c_cell_type_tx));
    [y] = ismember(upper(ctypes), upper(string(unique(T.CellType))));
    if ~any(y)
        errordlg('No primary cell type available in your data is supported by cellsubtype.xlsx.');
        return;
    end
    ctypelist = ctypes(y);

        if gui.i_isuifig(FigureHandle)
            [indx2, tf2] = gui.myListdlg(FigureHandle, cellstr(ctypelist), ...
                'Select Cell Type(s):');
        else
            [indx2, tf2] = listdlg('PromptString', 'Select Cell Type(s):', ...
                'SelectionMode', 'multiple', 'ListString', ...
                cellstr(ctypelist), 'ListSize', [220, 300]);
        end


        if tf2 ~= 1, return; end
        
        celltypetarget_list = ctypelist(indx2);

        answer = gui.myQuestdlg(FigureHandle, 'How to label cell type with subtype','Choose format', ...
            {'Type_{Subtype}','Type (Subtype)','Subtype'},'Type_{Subtype}');
        switch answer
            case 'Type_{Subtype}'
                formatid = 1;
            case 'Type (Subtype)'
                formatid = 2;
            case 'Subtype'
                formatid = 0;
        end

        for k = 1:length(celltypetarget_list)
            [sce] = sc_csubtypeanno(sce, celltypetarget_list(k), formatid);
        end
        if ~isa(src, "SingleCellExperiment")
            guidata(FigureHandle, sce);
        end
        requirerefresh = true;
end
