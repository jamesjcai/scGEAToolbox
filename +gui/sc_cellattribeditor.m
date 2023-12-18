function [sce] = sc_cellattribeditor(sce)

    T=table;
    if ~isempty(sce.list_cell_attributes)
        for k = 1:2:length(sce.list_cell_attributes)
            T = addvars(T, sce.list_cell_attributes{k+1}, ...
                'NewVariableNames', sce.list_cell_attributes{k});
        end
    end
    
    n=length(sce.list_cell_attributes);
    
    if ~isempty(sce.list_cell_attributes)
        for k = 1:2:length(sce.list_cell_attributes)
            attrlist = [attrlist; sce.list_cell_attributes{k}];
        end
    end
    
    [indx, tf] = listdlg('PromptString', {'Select a gene:', '', ''}, ...
        'SelectionMode', 'single', 'ListString', attrlist);
    
    x = inputdlg({'Attribute Name','Attribute Value'},...
                  'Attribute Editor', [1 50; 25 50]);

end
