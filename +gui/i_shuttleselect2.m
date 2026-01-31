function [idx] = i_shuttleselect2(items, preselected_items, parentfig)
    arguments
        items {mustBeValidItemList}
        preselected_items {mustBeValidItemList} = []
        parentfig = []
    end

    % Initial data setup
    if ~isempty(preselected_items)
        assert(all(ismember(preselected_items, items)));
        selItems = preselected_items(:);
        availItems = setdiff(items, selItems, 'stable');
    else
        selItems = {};
        availItems = items(:);
    end

    idx = [];
    hFig = uifigure('Name','Item Selection (Optimized)', ...
        'Position',[650 290 600 450],... % Increased height for search bar
        'CloseRequestFcn', @(src,evt)closeFcn(src),...
        'Visible','off');
    
    try gui.i_movegui2parent(hFig, parentfig); catch; end

    % --- Search Bar ---
    uilabel(hFig, 'Text', 'Search:', 'Position', [50 415 50 20]);
    efSearch = uieditfield(hFig, 'text', 'Position', [100 415 150 25], ...
        'Placeholder', 'Filter list...', ...
        'ValueChangedFcn', @(src, event) filterAvailable(src.Value));

    % --- Left Table: Available ---
    uilabel(hFig,'Text','Available Items','Position',[50 385 120 20]);
    tAvail = uitable(hFig, 'Position', [50 100 200 280]);
    tAvail.Data = table(availItems, 'VariableNames', {'Items'});
    setupTableLook(tAvail);

    % --- Right Table: Selected ---
    uilabel(hFig,'Text','Selected Items','Position',[350 385 120 20]);
    tSel = uitable(hFig, 'Position', [350 100 200 280]);
    tSel.Data = table(selItems, 'VariableNames', {'Items'});
    setupTableLook(tSel);

    % --- Buttons ---
    uibutton(hFig,'Text','>', 'Position',[270 300 60 30], 'ButtonPushedFcn',@(btn,evt) moveTableItems(tAvail, tSel));
    uibutton(hFig,'Text','>>', 'Position',[270 260 60 30], 'ButtonPushedFcn',@(btn,evt) moveAllTable(tAvail, tSel));
    uibutton(hFig,'Text','<', 'Position',[270 220 60 30], 'ButtonPushedFcn',@(btn,evt) moveTableItems(tSel, tAvail));
    uibutton(hFig,'Text','<<', 'Position',[270 180 60 30], 'ButtonPushedFcn',@(btn,evt) moveAllTable(tSel, tAvail));
    
    uibutton(hFig,'Text','Done', 'Position', [250 35 100 35], ...
        'ButtonPushedFcn', @(btn,evt) okTableFcn(hFig, tSel, items));   

    hFig.Visible = "on";
    uiwait(hFig);

    % --- Nested Callbacks ---
    function filterAvailable(query)
        % This filters based on the current 'availItems' state
        if isempty(query)
            tAvail.Data = table(availItems, 'VariableNames', {'Items'});
        else
            mask = contains(availItems, query, 'IgnoreCase', true);
            tAvail.Data = table(availItems(mask), 'VariableNames', {'Items'});
        end
    end

    function moveTableItems(src, dst)
        if isempty(src.Selection), return; end
        selRows = unique(src.Selection(:, 1));
        movedNames = src.Data.Items(selRows);
        
        % Update the master logical tracking
        if src == tAvail
            % Moving Available -> Selected
            selItems = [selItems; movedNames];
            availItems = setdiff(availItems, movedNames, 'stable');
        else
            % Moving Selected -> Available
            availItems = [availItems; movedNames];
            selItems = setdiff(selItems, movedNames, 'stable');
        end
        
        % Refresh UI
        filterAvailable(efSearch.Value);
        tSel.Data = table(selItems, 'VariableNames', {'Items'});
        src.Selection = [];
    end

    function moveAllTable(src, dst)
        if src == tAvail
            % Move only the CURRENTLY FILTERED items or literally everything?
            % Usually, 'Move All' implies moving the filtered subset. 
            % Here, we move all filtered items to selected:
            toMove = tAvail.Data.Items;
            selItems = [selItems; toMove];
            availItems = setdiff(availItems, toMove, 'stable');
        else
            toMove = tSel.Data.Items;
            availItems = [availItems; toMove];
            selItems = {};
        end
        filterAvailable(efSearch.Value);
        tSel.Data = table(selItems, 'VariableNames', {'Items'});
    end

    function okTableFcn(fig, t, originalItems)
        [~, idx] = ismember(selItems, originalItems);
        idx = idx(:);
        uiresume(fig);
        delete(fig);
    end
end

% --- Helper Functions ---
function setupTableLook(t)
    t.ColumnName = {}; t.RowName = {}; t.ColumnWidth = {'1x'};
    t.ColumnEditable = false; t.SelectionType = 'row';
end

function closeFcn(fig)
    uiresume(fig);
    delete(fig);
end

function mustBeValidItemList(x)
    if ~isempty(x) && ~(iscellstr(x) || isstring(x) || iscategorical(x) && isvector(x))
        error('Items must be a 1-D cell array of strings or similar.');
    end
end