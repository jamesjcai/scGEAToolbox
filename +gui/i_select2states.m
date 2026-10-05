function [thisc1, clabel1, thisc2, clabel2] = i_select2states(sce, ...
    allowsingle, parentfig)

if nargin<3, parentfig=[]; end
if nargin<2, allowsingle = false; end
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end

thisc1 = [];
clabel1 = '';
thisc2 = [];
clabel2 = '';


baselistitems = {'Current Class (C)'};
i_additem(sce.c_cluster_id, 'Cluster ID');
i_additem(sce.c_cell_cycle_tx, 'Cell Cycle Phase');
i_additem(sce.c_cell_type_tx, 'Cell Type');
i_additem(sce.c_batch_id, 'Batch ID');
i_additem(full(sum(sce.X))', 'Library Size');

function i_additem(itemv, itemn)
        if ~isempty(itemv) && length(unique(itemv)) >= 1
            baselistitems = [baselistitems, itemn];
        end
    end


listitems = [baselistitems, ...
        sce.list_cell_attributes(1:2:end)];
nx = length(baselistitems);

n = length(listitems);
if n < 2
       gui.myWarndlg(parentfig, ['This function requires at least two ', ...
                'grouping variables (e.g., BATCH_ID, ', ...
                'CLUSTER_ID, or CELL_TYPE_TXT).']);
        return;
    end


preftagname ='selected2states';
defaultindx = getpref('scgeatoolbox', preftagname, [n-1, n]);
if any(defaultindx > n) || any(defaultindx < 1), defaultindx = [n-1, n]; end

if gui.i_isuifig(parentfig)
        [indx2, tf2] = gui.myListdlg(parentfig, listitems, ...
            'Select cell state/grouping variable:', ...
            listitems(defaultindx), true);
    else
        [indx2, tf2] = listdlg('PromptString', ...
            {'Select cell state/grouping variable:'}, ...
            'SelectionMode', 'multiple', ...
            'ListString', listitems, ...
            'InitialValue', defaultindx, 'ListSize', [220, 300]);
    end

if tf2 == 1
        if length(indx2) ~= 2
            % Exactly one, when one is allowed. Three or more used to be
            % cut to the first without a word.
            if allowsingle && isscalar(indx2)
                [thisc1, clabel1] = i_getidx(indx2(1));
            elseif allowsingle
                gui.myWarndlg(parentfig, ...
                    'Please select one or two variables.','');
                return;
            else
                gui.myWarndlg(parentfig, ...
                    'Please select 2 grouping variables.','');
                return;
            end
        else
            [thisc1, clabel1] = i_getidx(indx2(1));
            [thisc2, clabel2] = i_getidx(indx2(2));
        end
        setpref('scgeatoolbox', preftagname, indx2);
    end


function [thisc, clabel] = i_getidx(indx)
        clabel = listitems{indx};
        switch clabel
            % Full double columns: SCE.X is single sparse from R2025a on,
            % and STRING() cannot convert a sparse array. See the note in
            % GUI.I_SELECT1STATE.
            case 'Library Size'
                thisc = full(double(sum(sce.X, 1)))';
            case 'Mt-reads Ratio'
                i = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
                lbsz = full(double(sum(sce.X, 1)));
                lbsz_mt = full(double(sum(sce.X(i, :), 1)));
                thisc = (lbsz_mt ./ lbsz).';
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
            otherwise % other properties
                nx = length(baselistitems);
                clabel = sce.list_cell_attributes{2*(indx - nx)-1};
                thisc = sce.list_cell_attributes{2*(indx - nx)};
        end
    end
end
