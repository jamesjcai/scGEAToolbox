function [thisc1, clabel1, thisc2, clabel2] = i_select2class(sce, ...
    allowsingle)

if nargin<2, allowsingle = false; end


thisc1 = [];
clabel1 = '';
thisc2 = [];
clabel2 = '';

listitems = {'Current Class (C)'};
i_additem(sce.c_cluster_id, 'Cluster ID');
i_additem(sce.c_cell_cycle_tx, 'Cell Cycle Phase');
i_additem(sce.c_cell_type_tx, 'Cell Type');
i_additem(sce.c_batch_id, 'Batch ID');

function i_additem(itemv, itemn)
    if ~isempty(itemv) && length(unique(itemv)) > 1
        listitems = [listitems, itemn];
    end
end


    n = length(listitems);
    if n < 2
        errordlg('Need at least two grouping variables (e.g., BATCH_ID, CLUSTER_ID, or CELL_TYPE_TXT)');
        return;
    end

    % listitems={'Current Class (C)','Cluster ID','Batch ID',...
    %            'Cell Type','Cell Cycle Phase'};
    [indx2, tf2] = listdlg('PromptString', ...
        {'Select two grouping varibles:'}, ...
        'SelectionMode', 'multiple', ...
        'ListString', listitems, ...
        'InitialValue', [n - 1, n], 'ListSize', [220, 300]);
    if tf2 == 1
        if length(indx2) ~= 2
            if allowsingle
                [thisc1, clabel1] = i_getidx(indx2(1));
            else
                warndlg('Please select 2 grouping variables.','');
                return;                
            end
        else
            [thisc1, clabel1] = i_getidx(indx2(1));
            [thisc2, clabel2] = i_getidx(indx2(2));
        end
    end

    % if isempty(thisc)
    %     errordlg('Undefined');
    %     return;
    % end
    % if numel(unique(thisc))==1
    %     warndlg("Cannot compare with an unique group");
    %     return;
    % end


    function [thisc, clabel] = i_getidx(indx)
        clabel = listitems{indx};
        switch clabel
            case 'Current Class (C)'
                thisc = sce.c;
            case 'Cluster ID' % cluster id
                thisc = sce.c_cluster_id;
            case 'Batch ID' % batch id
                thisc = sce.c_batch_id;
            case 'Cell Type' % cell type
                thisc = sce.c_cell_type_tx;
            case 'Cell Cycle Phase' % cell cycle
                thisc = sce.c_cell_cycle_tx;
        end
    end

end
