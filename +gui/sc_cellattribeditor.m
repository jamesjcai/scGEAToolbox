function [sce, needupdate] = sc_cellattribeditor(sce, addnew)
if nargin<2, addnew = false; end

    needupdate = false;

if ~addnew    
    baselistitems = {'Cell Cycle Phase', ...
        'Cell Type', 'Cluster ID', ...
        'Batch ID', 'Cell ID'};
    listitems = [baselistitems, ...
        sce.list_cell_attributes(1:2:end)];
    listitems = listitems(~cellfun(@isempty, listitems));

    [indx2, tf2] = listdlg('PromptString', ...
        {'Select Cell Attribute:'}, ...
        'SelectionMode', 'single', 'ListString', listitems,'ListSize',[200 300]);
    
    if tf2 == 1
        clable = listitems{indx2};
        switch clable
            case 'Cluster ID'   % cluster id
                thisc = sce.c_cluster_id;
            case 'Batch ID'       % batch id
                thisc = sce.c_batch_id;
            case 'Cell Type'  % cell type
                thisc = sce.c_cell_type_tx;
            case 'Cell Cycle Phase' % cell cycle
                thisc = sce.c_cell_cycle_tx;
            case 'Cell ID'
                thisc = sce.c_cell_id;
            otherwise
                [y,idx]=ismember(clable, sce.list_cell_attributes(1:2:end));
                if y
                    thisc = sce.list_cell_attributes{idx+1};
                end
        end
    else
        return;
    end
    
    %fw = gui.gui_waitbar;
    x = inputdlg(sprintf('Attribute Name: %s\n%s',clable, 'Attribute Value:'), ...
                      'Attribute Editor', [25 50], {char(string(thisc))});
    %gui.gui_waitbar(fw);

else  
    x = inputdlg({'Attribute Name','Attribute Value'},...
                  'Attribute Editor', [1 50; 25 50]);
    % {'new_attrib', char(string([1:sce.NumCells]'))});
end

    if isempty(x), return; end


    answer = questdlg('What is the input data type?', ...
	    'Data Type', ...
	    'String','Numeric','Cancel','String');
    switch answer
        case 'String'
            if addnew
                clable = strtrim(string(x{1}));
                newthisc = strtrim(string(x{2}));
            else
                newthisc = strtrim(string(x{1}));
            end
        case 'Numeric'
            if addnew
                clable = strtrim(string(x{1}));
                newthisc = str2double(string(x{2}));
            else
                newthisc = str2double(string(x{1}));
            end
        case 'Cancel'
            return;
    end


    if addnew
        if size(newthisc,1)~=sce.NumCells
           errordlg('Attribute length is not equal to the number of cells.');
           return;
        end
    else
        if ~isequal(size(newthisc), size(thisc))
           errordlg('Attribute length changed.');
           return;
        end
    end

    if addnew
        clable = matlab.lang.makeValidName(clable);
        sce.list_cell_attributes = [sce.list_cell_attributes, ...
                                    {clable, newthisc(:)}];
        waitfor(helpdlg('Cell Attribute Added.',''));
        needupdate = true;
    else
        switch clable
            case 'Cluster ID'   % cluster id
                sce.c_cluster_id = newthisc;
            case 'Batch ID'       % batch id
                sce.c_batch_id = newthisc;
            case 'Cell Type'  % cell type
                sce.c_cell_type_tx = newthisc;
            case 'Cell Cycle Phase' % cell cycle
                sce.c_cell_cycle_tx = newthisc;
            case 'Cell ID'
                sce.c_cell_id = newthisc;
            otherwise
                [y,idx]=ismember(clable, sce.list_cell_attributes(1:2:end));
                if y, sce.list_cell_attributes{idx+1} = newthisc; end
        end
        waitfor(helpdlg('Cell Attribute Changed.',''));
        needupdate = true;
    end    
end