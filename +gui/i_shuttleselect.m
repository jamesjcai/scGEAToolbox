function [idx] = i_shuttleselect(items, preselected_items, parentfig)
    % Optimized for large datasets (10k to 1M+ items)
    if ~isempty(preselected_items)
        selItems = preselected_items(:);
        availItems = setdiff(items, selItems, 'stable');
    else
        selItems = {};
        availItems = items(:);
    end

    idx = [];
    hFig = uifigure('Name','Item Selection (High Capacity)', ...
        'Position',[650 290 600 450], 'Visible','off');
    
    % Internal helper to configure table appearance
    function applyStyle(t)
        t.ColumnName = {}; 
        t.RowName = {}; 
        t.ColumnWidth = {'1x'};
        t.ColumnEditable = false; 
        t.SelectionType = 'row';
    end

    try gui.i_movegui2parent(hFig, parentfig); catch; end

    % --- Search Bar Row ---
    uilabel(hFig, 'Text', 'Search:', 'Position', [50 415 50 20]);
    efSearch = uieditfield(hFig, 'text', 'Position', [100 415 150 25], ...
        'Placeholder', 'Filter list...', ...
        'ValueChangedFcn', @(src, event) filterTable());
    
    uibutton(hFig, 'Text', 'X', 'Position', [255 415 25 25], ...
        'Tooltip', 'Clear Search', ...
        'ButtonPushedFcn', @(btn, event) clearSearch());

    % --- Status Labels (Counters) ---
    lblAvailCount = uilabel(hFig, 'Text', '', 'Position', [50 80 200 20], ...
        'FontAngle', 'italic', 'FontColor', [0.4 0.4 0.4]);
    lblSelCount = uilabel(hFig, 'Text', '', 'Position', [350 80 200 20], ...
        'FontAngle', 'italic', 'FontColor', [0.4 0.4 0.4]);

    % --- Tables ---
    uilabel(hFig,'Text','Available Items','Position',[50 385 120 20], 'FontWeight', 'bold');
    tAvail = uitable(hFig, 'Position', [50 100 200 280]);
    
    uilabel(hFig,'Text','Selected Items','Position',[350 385 120 20], 'FontWeight', 'bold');
    tSel = uitable(hFig, 'Position', [350 100 200 280]);
    
    applyStyle(tAvail);
    applyStyle(tSel);
    
    % Initial load
    refreshUI();

    % --- Shuttle Buttons ---
    uibutton(hFig,'Text','>', 'Position',[270 300 60 30], 'ButtonPushedFcn',@(btn,evt) moveTable(tAvail, tSel, true));
    uibutton(hFig,'Text','>>', 'Position',[270 260 60 30], 'ButtonPushedFcn',@(btn,evt) moveAllTable(tAvail, tSel, true));
    uibutton(hFig,'Text','<', 'Position',[270 220 60 30], 'ButtonPushedFcn',@(btn,evt) moveTable(tSel, tAvail, false));
    uibutton(hFig,'Text','<<', 'Position',[270 180 60 30], 'ButtonPushedFcn',@(btn,evt) moveAllTable(tSel, tAvail, false));
    
    uibutton(hFig,'Text','Done', 'Position', [250 35 100 35], 'FontWeight', 'bold', ...
        'ButtonPushedFcn', @(btn,evt) uiresume(hFig));    

    hFig.Visible = "on";
    uiwait(hFig);

    if isvalid(hFig)
        [~, idx] = ismember(selItems, items);
        idx = idx(:);
        delete(hFig);
    end

    % --- Nested Callbacks ---
    function clearSearch()
        efSearch.Value = '';
        refreshUI();
    end

    function filterTable()
        query = efSearch.Value;
        if isempty(query)
            tAvail.Data = table(availItems, 'VariableNames', {'Items'});
        else
            mask = contains(availItems, query, 'IgnoreCase', true);
            tAvail.Data = table(availItems(mask), 'VariableNames', {'Items'});
        end
        updateCounters();
    end

    function moveTable(src, dst, isAdding)
        if isempty(src.Selection), return; end
        rows = unique(src.Selection(:, 1));
        names = src.Data.Items(rows);
        if isAdding
            selItems = [selItems; names];
            availItems = setdiff(availItems, names, 'stable');
        else
            availItems = [availItems; names];
            selItems = setdiff(selItems, names, 'stable');
        end
        refreshUI();
    end

    function moveAllTable(src, dst, isAdding)
        names = src.Data.Items; 
        if isempty(names), return; end
        if isAdding
            selItems = [selItems; names];
            availItems = setdiff(availItems, names, 'stable');
        else
            availItems = [availItems; names];
            selItems = setdiff(selItems, names, 'stable');
        end
        refreshUI();
    end

    function updateCounters()
        numVisible = height(tAvail.Data);
        numTotalAvail = numel(availItems);
        if isempty(efSearch.Value)
            lblAvailCount.Text = sprintf('%d items available', numTotalAvail);
        else
            lblAvailCount.Text = sprintf('Found %d of %d items', numVisible, numTotalAvail);
        end
        lblSelCount.Text = sprintf('%d items selected', numel(selItems));
    end

    function refreshUI()
        filterTable();
        tSel.Data = table(selItems, 'VariableNames', {'Items'});
        tAvail.Selection = []; tSel.Selection = [];
        updateCounters();
    end
end