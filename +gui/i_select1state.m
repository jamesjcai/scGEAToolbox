function [thisc, clabel, listitems, newpickclabel] = i_select1state(sce, ...
    nobaseline, nocustome, noattrib, continuousonly, parentfig)

if nargin < 2, nobaseline = false; end
if nargin < 3, nocustome = false; end
if nargin < 4, noattrib = true; end
if nargin < 5, continuousonly = false; end
if nargin < 6, parentfig = []; end

if continuousonly
end

thisc = [];
clabel = '';
newpickclabel = '';

if nobaseline
    listitems = sce.list_cell_attributes(1:2:end);
else
    baselistitems = {'Library Size', ...
        'Number of Detected Genes', ...
        'Mt-reads Ratio', ...
        '-------------------', ...
        'Cell Cycle Phase', ...
        'Cell Type', ...
        'Cluster ID', ...
        'Batch ID', ...
        '-------------------'};
    listitems = [baselistitems, ...
        sce.list_cell_attributes(1:2:end)];
end
listitems = listitems(~cellfun(@isempty, listitems));

if ~nocustome
    a = evalin('base', 'whos');
    b = struct2cell(a);
    v = false(length(a), 1);
    for k = 1:length(a)
        if max(a(k).size) == sce.NumCells && min(a(k).size) == 1 && ...
                ~strcmp(a(k).class, 'table')
            v(k) = true;
        end
    end
    if any(v)
        a = a(v);
        b = b(:, v);
        listitems = [listitems, 'Workspace Variable...'];
    end
end

if ~noattrib
    % if istable(sce.table_attributes)
    %     if size(sce.table_attributes, 1) == sce.NumCells
    %         listitems = [listitems, 'SCE Attribute Table...'];
    %     end
    % end
end

%if ~(ismcc || isdeployed)
listitems = [listitems, 'Load Variable from File...'];
%end

if isempty(listitems), return; end


% listitems={'Current Class (C)','Cluster ID','Batch ID',...
%            'Cell Type','Cell Cycle Phase'};


 if gui.i_isuifig(parentfig)
    [indx2, tf2] = gui.myListdlg(parentfig, listitems, 'Select state/grouping variable:');
 else
    [indx2, tf2] = listdlg('PromptString', ...
        {'Select state/grouping variable:'}, ...
        'SelectionMode', 'single', 'ListString', listitems, ...
        'ListSize', [220, 300]);
 end
if tf2 == 1
    clabel = listitems{indx2};
    switch clabel
        case 'Library Size'
            thisc = sum(sce.X).';
        case 'Number of Detected Genes'
            thisc = sum(sce.X > 0, 1);
        case 'Mt-reads Ratio'
            i = startsWith(sce.g, 'mt-', 'IgnoreCase', true);
            lbsz = sum(sce.X, 1);
            lbsz_mt = sum(sce.X(i, :), 1);
            thisc = (lbsz_mt ./ lbsz).';
        case 'Cluster ID' % cluster id
            thisc = sce.c_cluster_id;
        case 'Batch ID' % batch id
            thisc = sce.c_batch_id;
        case 'Cell Type' % cell type
            thisc = sce.c_cell_type_tx;
        case 'Cell Cycle Phase' % cell cycle
            needestimate = false;
            if isempty(sce.c_cell_cycle_tx) || all(strcmp(unique(sce.c_cell_cycle_tx), "undetermined"))
                needestimate = true;
            else
                answer1 = gui.myQuestdlg(parentfig, 'Use existing cell cycle estimation or re-compute new estimation?', ...
                    '', {'Use existing', 'Re-compute', 'Cancel'}, 'Use existing');
                switch answer1
                    case 'Re-compute'
                        needestimate = true;
                    case 'Cancel'
                        return;
                end
            end
            if needestimate
                fw = gui.myWaitbar(parentfig);
                sce = sce.estimatecellcycle(true, 1);
                gui.myWaitbar(parentfig, fw);
            end
            thisc = sce.c_cell_cycle_tx;
        case '-------------------'
            thisc = [];
            clabel = '';
            listitems = [];
            newpickclabel = [];
            return;
        case 'Workspace Variable...'
            [thisc, newpickclabel] = i_pickvariable;
        case 'SCE Attribute Table...'
            [thisc, newpickclabel] = i_pickattribute;
        case 'Load Variable from File...'
            thisc = [];
            clabel = '';
            listitems = [];
            newpickclabel = [];

            [fname, pathname] = uigetfile( ...
                {'*.txt', 'Variable Data Files (*.txt)'; ...
                '*.*', 'All Files (*.*)'}, ...
                'Pick a Variable Data File');
            if isvalid(parentfig) && isa(parentfig, 'matlab.ui.Figure'), figure(parentfig); end
            if isequal(fname, 0), return; end
            txtfile = fullfile(pathname, fname);
            try
                T = readtable(txtfile);
                thisc = T.(T.Properties.VariableNames{1});
                assert(length(thisc) == sce.NumCells, 'Invalid Variable Data File.')
                newpickclabel = 'ExternalVariable';
            catch ME
                thisc = [];
                gui.myErrordlg(parentfig, ME.message, ME.identifier);
                return;
            end
        otherwise % other properties

            [~, idx] = ismember(clabel, ...
                string(sce.list_cell_attributes(1:2:end)));
            thisc = sce.list_cell_attributes{idx*2};

            %nx=length(baselistitems);
            %clabel = sce.list_cell_attributes{2 * (indx2 - nx) - 1};
            %thisc = sce.list_cell_attributes{2 * (indx2 - nx)};
    end
end


    function [c, x] = i_pickvariable
        c = [];
        x = '';
        %     a=evalin('base','whos');
        %     b=struct2cell(a);
        %     v=false(length(a),1);
        %     for k=1:length(a)
        %         if max(a(k).size)==sce.NumCells && min(a(k).size)==1
        %             v(k)=true;
        %         end
        %     end
        %     if any(v)
        %valididx=ismember(b(4,:),'double');
        %a=a(valididx);
        
        if gui.i_isuifig(parentfig)
            [indx, tf] = gui.myListdlg(parentfig, b(1, :), 'Select workspace variable:');
        else
            [indx, tf] = listdlg('PromptString', {'Select workspace variable:'}, ...
                'liststring', b(1, :), 'SelectionMode', 'single', 'ListSize', [220, 300]);
        end

        if tf == 1
            c = evalin('base', a(indx).name);
            x = a(indx).name;
        end
        %    end
end

        function [c, x] = i_pickattribute
            c = [];
            x = '';
            T = sce.table_attributes;
            att = sce.table_attributes.Properties.VariableNames;
            
        if gui.i_isuifig(parentfig)
            [indx, tf] = gui.myListdlg(parentfig, att, ...
                'Select a SCE attribute variable:');
        else
            [indx, tf] = listdlg('PromptString', {'Select a SCE attribute variable:'}, ...
                'liststring', att, 'SelectionMode', 'single', 'ListSize', [220, 300]);
        end

            if tf == 1
                x = string(att(indx));
                c = T.(x);
                if iscell(c)
                    c = string(c);
                end
            end
    end


end