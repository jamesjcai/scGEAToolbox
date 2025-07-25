function [sce, needupdate] = sc_cellattribeditor(sce, addnew, parentfig)
if nargin<3, parentfig = []; end
if nargin<2, addnew = false; end

    needupdate = false;


if ~addnew    % edit
    baselistitems = {'Cell Cycle Phase', ...
        'Cell Type', 'Cluster ID', ...
        'Batch ID', 'Cell ID'};
    listitems = [baselistitems, ...
        sce.list_cell_attributes(1:2:end)];
    listitems = listitems(~cellfun(@isempty, listitems));

        if gui.i_isuifig(parentfig)
            [indx2, tf2] = gui.myListdlg(parentfig, listitems, 'Select Cell Attribute:');
        else
            [indx2, tf2] = listdlg('PromptString', ...
                {'Select Cell Attribute:'}, ...
                'SelectionMode', 'single', 'ListString', listitems, ...
                'ListSize', [220, 300]);
        end
    
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
                [y,idx] = ismember(clabel, sce.list_cell_attributes(1:2:end));
                if y
                    thisc = sce.list_cell_attributes{2*idx};
                end
        end
    else
        return;
    end
    
  
    
%    if ~strcmp('Yes', gui.myQuestdlg(parentfig, ...
%            'It may take a while to load values. Continue?'))
%        return; 
%    end
    tic;
    if gui.i_isuifig(parentfig)
        %        x = gui.myInputdlg({sprintf('Attribute Name: %s\n%s',clabel, 'Attribute Values:')}, ...
        %                          'Attribute Editor', {char(string(thisc))}, parentfig);
        
        % assignin("base","thisc",thisc);
        % assignin("base","clabel",clabel);

        x = gui.myTextareadlg(parentfig, {'Attribute Name', 'Attribute Values'},...
                      'Attribute Editor', {clabel, string(thisc)}, [false, true]);
        % assignin("base","x",x);
        if ~isempty(x)
            x(1)=[];
        end
    else
        x = inputdlg(sprintf('Attribute Name: %s\n%s',clabel, 'Attribute Values:'), ...
                          'Attribute Editor', [15 80], {char(string(thisc))});
        
    end
    toc;


else    % add new

    if gui.i_isuifig(parentfig)
        % x = gui.myInputdlg({'Attribute Name','Attribute Values'},...
        %              'Attribute Editor', {''}, parentfig); % Assuming default is empty cell
        % x = gui.myTextareadlg(parentfig, '', 'Attribute Name'); % Assuming default is empty cell
        x = gui.myTextareadlg(parentfig, {'Attribute Name','Attribute Values'},...
                      'Attribute Editor', {'New_Attribute', ("Value_"+string(1:sce.NumCells))'});
    else
        x = inputdlg({'Attribute Name','Attribute Values'},...
                      'Attribute Editor', [1 80; 15 80]);
    end
    % {'new_attrib', char(string([1:sce.NumCells]'))});
end

if isempty(x), return; end

if addnew
    if isempty(x{1})
        gui.myWarndlg(parentfig, 'Attribute Name cannot be empty.');
        return; 
    end
    if isempty(x{2})
        gui.myWarndlg(parentfig, 'Attribute Values cannot be empty.');
        return; 
    end
else
    if isempty(x{1})       % when add new - x{1} is the values
        gui.myWarndlg(parentfig, 'Attribute Values cannot be empty.');
        return;
    end
end

    answer = gui.myQuestdlg(parentfig, 'What is the data type of attribute values?', ...
	    'Data Type', ...
	    {'String','Numeric','Cancel'},'String');
    switch answer
        case 'String'
            if addnew
                clabel = strtrim(x{1});
                newthisc = strtrim(string(trimBottomEmpty(x{2})));
            else
                newthisc = strtrim(string(trimBottomEmpty(x{1})));    % when add new - x{1} is the values
            end
        case 'Numeric'
            if addnew
                clabel = strtrim(x{1});
                newthisc = str2double(string(trimBottomEmpty(x{2})));
            else
                newthisc = str2double(string(trimBottomEmpty(x{1})));   % when add new - x{1} is the values
            end
        case 'Cancel'
            return;
    end


    if addnew
        if size(newthisc,1) ~= sce.NumCells
           gui.myWarndlg(parentfig, ...
               'Attribute length is not equal to the number of cells.');
           return;
        end
    else
        if ~isequal(size(newthisc), size(thisc))
           gui.myWarndlg(parentfig, 'Attribute length changed.');
           return;
        end
    end

    if addnew
        clabel = matlab.lang.makeValidName(clabel);        
        existinglabels = sce.list_cell_attributes(1:2:end);
        if ismember(clabel, existinglabels)
            gui.myWarndlg(parentfig, 'Cell Attribute Name Existing.');
            return;
        end
        sce.list_cell_attributes = [sce.list_cell_attributes, ...
                                    {clabel, newthisc(:)}];
        gui.myHelpdlg(parentfig, 'Cell Attribute Added.');
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
        gui.myHelpdlg(parentfig, 'Cell Attribute Changed.');
        needupdate = true;
    end    
end




function trimmed_text = trimBottomEmpty(input_text)
    % TRIMBOTTOMEMPTY Trims whitespace from each line and removes empty lines from bottom
    %
    % SYNTAX:
    %   trimmed_text = trimBottomEmpty(input_text)
    %
    % INPUT:
    %   input_text - char array or string, can be multiline
    %
    % OUTPUT:
    %   trimmed_text - char array with whitespace trimmed from each line
    %                  and empty lines removed from bottom only
    %
    % EXAMPLE:
    %   text = ['111'; '222'; '333'; '   '];
    %   result = trimBottomEmpty(text);
    
    lines = cellstr(input_text);     % Convert to cell array
    lines = strtrim(lines);          % Trim each line
    
    % Find the last non-empty line
    last_non_empty = find(~cellfun('isempty', lines), 1, 'last');
    
    if ~isempty(last_non_empty)
        lines = lines(1:last_non_empty);  % Keep only up to last non-empty line
    else
        lines = {};  % All lines were empty
    end
    
    trimmed_text = char(lines);      % Convert back to char array
end


