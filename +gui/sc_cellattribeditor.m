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
        'SelectionMode', 'single', 'ListString', listitems, ...
        'ListSize', [220, 300]);
    
    if tf2 == 1
        clabel = listitems{indx2};
        switch clabel
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
                [y,idx]=ismember(clabel, sce.list_cell_attributes(1:2:end));
                if y
                    thisc = sce.list_cell_attributes{idx+1};
                end
        end
    else
        return;
    end
    
    %fw = gui.gui_waitbar;
    uiwait(helpdlg('It may take a while to load values. Click OK and wait.',''));
    x = inputdlg(sprintf('Attribute Name: %s\n%s',clabel, 'Attribute Values:'), ...
                      'Attribute Editor', [15 80], {char(string(thisc))});
    %gui.gui_waitbar(fw);

else  
    x = inputdlg({'Attribute Name','Attribute Values'},...
                  'Attribute Editor', [1 80; 15 80]);
    % {'new_attrib', char(string([1:sce.NumCells]'))});
end

if isempty(x), return; end

if addnew
    if isempty(x{1})
        warndlg('Attribute Name cannot be empty.','','modal');
        return; 
    end
    if isempty(x{2})
        warndlg('Attribute Values cannot be empty.','','modal');
        return; 
    end
else
    if isempty(x{1})       % when add new - x{1} is the values
        warndlg('Attribute Values cannot be empty.','','modal');
        return;
    end
end

    answer = questdlg('What is the data type of attribute values?', ...
	    'Data Type', ...
	    'String','Numeric','Cancel','String');
    switch answer
        case 'String'
            if addnew
                clabel = strtrim(x{1});
                newthisc = strtrim(string(x{2}));
            else
                newthisc = strtrim(string(x{1}));    % when add new - x{1} is the values
            end
        case 'Numeric'
            if addnew
                clabel = strtrim(x{1});
                newthisc = str2double(string(x{2}));
            else
                newthisc = str2double(string(x{1}));   % when add new - x{1} is the values
            end
        case 'Cancel'
            return;
    end


    if addnew
        if size(newthisc,1)~=sce.NumCells
           warndlg('Attribute length is not equal to the number of cells.','','modal');
           return;
        end
    else
        if ~isequal(size(newthisc), size(thisc))
           warndlg('Attribute length changed.','','modal');
           return;
        end
    end

    if addnew
        clabel = matlab.lang.makeValidName(clabel);        
        existinglabels = sce.list_cell_attributes(1:2:end);
        if ismember(clabel, existinglabels)
            warndlg('Cell Attribute Name Existing.','','modal');
            return;
        end
        sce.list_cell_attributes = [sce.list_cell_attributes, ...
                                    {clabel, newthisc(:)}];
        waitfor(helpdlg('Cell Attribute Added.',''));
        needupdate = true;
    else
        switch clabel
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
                [y,idx]=ismember(clabel, sce.list_cell_attributes(1:2:end));
                if y, sce.list_cell_attributes{idx+1} = newthisc; end
        end
        waitfor(helpdlg('Cell Attribute Changed.',''));
        needupdate = true;
    end    
end
