function [indx, tf] = myTabledlg(parentfig, options, Title, ...
        prefersel, allowmulti)
    if nargin < 5, allowmulti = true; end
    if nargin < 4, prefersel = []; end
    
    parentPos = parentfig.Position;
    parentCenter = [parentPos(1) + parentPos(3)/2, parentPos(2) + parentPos(4)/2];
    dlgSize = [300, 450]; 
    dlgPos = [parentCenter(1) - dlgSize(1)/2, parentCenter(2) - dlgSize(2)/2, dlgSize];

    d = uifigure('Name', Title, 'Position', dlgPos, ...
        'WindowStyle', 'normal', 'Visible','off');
    
    % --- UITable Setup (Replacing Listbox) ---
    % uitable doesn't have a 'Value' property for strings, so we manage indices
    ut = uitable(d, 'Position', [20 60 260 370]);
    ut.Data = table(options(:), 'VariableNames', {'Items'});
    ut.ColumnName = {};
    ut.RowName = {};
    ut.ColumnWidth = {'1x'};
    ut.ColumnEditable = false;
    
    if allowmulti
        ut.SelectionType = 'row'; % Allows shift/ctrl click
    else
        ut.SelectionType = 'row'; 
        % Note: uitable doesn't have a strict 'single' mode like listbox, 
        % but we can enforce it in the SelectionChangedFcn if absolutely necessary.
    end

    % Handle pre-selection
    if ~isempty(prefersel)
        [~, pre_idx] = ismember(prefersel, options);
        % In uitable, Selection is [row, col]. For whole row selection:
        ut.Selection = [pre_idx(:), ones(numel(pre_idx), 1)];
    end

    % Update KeyPress to work with table rows
    d.KeyPressFcn = @(src, event) jumpToFirstMatchTable(ut, options, event);

    % Buttons
    btnOK = uibutton(d, 'Text', 'OK', 'Position', [60 20 80 30], ...
        'ButtonPushedFcn', @(btn,event) uiresume(d));
    btnCancel = uibutton(d, 'Text', 'Cancel', 'Position', [160 20 80 30], ...
        'ButtonPushedFcn', @(btn,event) (close(d)));

    if ~isMATLABReleaseOlderThan('R2025a')
        try theme(d, parentfig.Theme.BaseColorStyle); catch; end
    end

    drawnow;
    d.Visible = 'on';
    d.WindowStyle = "modal";
    focus(ut); 

    uiwait(d);
    
    % --- Process Results ---
    if isvalid(ut)
        % Selection property returns N-by-2 matrix [row, col]
        rows = ut.Selection(:, 1);
        if isempty(rows)
            indx = [];
            tf = 0;
        else
            indx = unique(rows); % Ensure unique indices if user clicked weirdly
            tf = 1;
        end
        delete(d);
    else
        indx = [];
        tf = 0;
    end
end

function jumpToFirstMatchTable(ut, options, event)
    key = event.Character;
    if isempty(key) || ~ischar(key), return; end
    
    idx = find(startsWith(options, key, 'IgnoreCase', true), 1);
    if ~isempty(idx)
        % Move selection to the first match
        ut.Selection = [idx, 1];
        % Scroll to the row to ensure it's visible (Available in newer MATLAB versions)
        scroll(ut, 'row', idx);
    end
end