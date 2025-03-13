function [thisc1, clabel1, thisc2, clabel2] = i_select2class(sce, ...
    allowsingle, parentfig)

if nargin<3, parentfig=[]; end
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
        gui.myErrordlg(parentfig, 'Need at least two grouping variables (e.g., BATCH_ID, CLUSTER_ID, or CELL_TYPE_TXT)');
        return;
    end

    % listitems={'Current Class (C)','Cluster ID','Batch ID',...
    %            'Cell Type','Cell Cycle Phase'};

        if gui.i_isuifig(parentfig)
            [indx2, tf2] = gui.myListdlg(parentfig, listitems, 'Select two grouping varibles:');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select two grouping varibles:'}, ...
                'SelectionMode', 'multiple', ...
                'ListString', listitems, ...
                'InitialValue', [n - 1, n], 'ListSize', [220, 300]);
        end

    if tf2 == 1
        if length(indx2) ~= 2
            if allowsingle
                [thisc1, clabel1] = i_getidx(indx2(1));
            else
                gui.myWarndlg(parentfig, ...
                    'Please select 2 grouping variables.','');
                return;                
            end
        else
            [thisc1, clabel1] = i_getidx(indx2(1));
            [thisc2, clabel2] = i_getidx(indx2(2));
        end
    end

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
