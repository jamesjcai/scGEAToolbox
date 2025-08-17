function idx = i_selmultidialog(items, preselected_items, parentfig)
    % selectionDialog creates a dialog for selecting multiple items.
    %
    % Usage:
    %   idx = selectionDialog(items, preselected_items)
    %
    % Inputs:
    %   items              - cell array of strings (all available items)
    %   preselected_items  - cell array of strings (initial selected items)
    %
    % Output:
    %   idx - indices of selected items in `items`

%{
items = arrayfun(@(x) sprintf('Item %d',x), 1:10, 'UniformOutput', false);
pre = {'Item 3','Item 7'};

idx = gui.i_selmultidialog(items, pre)
idx2 = gui.i_selmultidialog(items, pre)

%}

    % arguments
    %     items cell
    %     preselected_items cell = {}
    %     parentfig
    % end
    % arguments
    %     items cell
    %     preselected_items cell = {}
    %     parentfig matlab.ui.Figure = []
    % end

arguments
    items {mustBeValidItemList}
    preselected_items {mustBeValidItemList} = []
    parentfig = []   % allow [] or matlab.ui.Figure    
end

    if ~isempty(preselected_items)
        assert(all(ismember(preselected_items, items)));
    end

    % Validate parentfig if non-empty
    if ~isempty(parentfig) && ~isa(parentfig,'matlab.ui.Figure')
        error('parentfig must be a matlab.ui.Figure or [].');
    end
    
    % if nargin < 2
    %     preselected_items = {};
    % end
    % if nargin < 3
    %     parentfig = [];
    % end

    if ~isempty(preselected_items)
        selItems = preselected_items(:)';
        availItems = setdiff(items, selItems, 'stable');
    else
        selItems = {};
        availItems = items;
    end

    idx = [];

    % UI figure
    hFig = uifigure('Name','Item Selection', ...
        'Position',[100 100 600 400],...
        'CloseRequestFcn',@(src,evt)closeFcn(src));
    gui.i_movegui2parent(hFig, parentfig);


    % Available list
    uilabel(hFig,'Text','Available Items','Position',[50 360 120 20]);
    lbAvail = uilistbox(hFig,...
        'Position',[50 100 200 260],...
        'Items',availItems,...
        'Multiselect','on');

    % Selected list
    % if isempty(selItems), selItems = {}; end
    uilabel(hFig,'Text','Selected Items','Position',[350 360 120 20]);
    lbSel = uilistbox(hFig,...
        'Position',[350 100 200 260],...
        'Items',selItems,...
        'Multiselect','on');

    % Buttons
    uibutton(hFig,'Text','>',...
        'Position',[270 280 60 30],...
        'ButtonPushedFcn',@(btn,event) moveItems(lbAvail,lbSel));
    uibutton(hFig,'Text','>>',...
        'Position',[270 240 60 30],...
        'ButtonPushedFcn',@(btn,event) moveAll(lbAvail,lbSel));
    uibutton(hFig,'Text','<',...
        'Position',[270 200 60 30],...
        'ButtonPushedFcn',@(btn,event) moveItems(lbSel,lbAvail));
    uibutton(hFig,'Text','<<',...
        'Position',[270 160 60 30],...
        'ButtonPushedFcn',@(btn,event) moveAll(lbSel,lbAvail));

    btnDone = uibutton(hFig,'Text','Done',...
        'Position',[250 40 100 30],...
        'ButtonPushedFcn',@(btn,event) okFcn(hFig,lbSel));    
    % Done button
    % btnDone = uibutton(hFig,'Text','Done',...
    %     'Position',[250 40 100 30],...
    %     'ButtonPushedFcn',@(btn,event) uiresume(hFig));
    % 'ButtonPushedFcn',@(btn,evt)okFcn(fig,lb)
    % Wait for user
    uiwait(hFig);

    function okFcn(fig,lb)
        % Store selection (as row) and close
        val = lb.Items;
        [~, idx] = ismember(val, items);
        idx=idx(:);
        uiresume(fig);
        delete(fig);
    end
    % Collect result
    %selItems = lbSel.Items;
    %[~, idx] = ismember(selItems, items);
    %idx = idx(:);

    % Close the figure
    %delete(hFig);
    if isvalid(hFig), delete(hFig); end
end

% --- Helper functions ---

    function closeFcn(fig)
        uiresume(fig);
        delete(fig);
    end

function moveItems(srcList, dstList)
    items = srcList.Items;
    sel   = srcList.Value;

    % Remove from source
    srcList.Items(ismember(items,sel)) = [];

    % Add to destination
    dstList.Items = [dstList.Items(:); sel(:)];
end

function moveAll(srcList, dstList)
    dstList.Items = [dstList.Items(:); srcList.Items(:)];
    srcList.Items = {};
end

%{
function i_selmultidialog(genelist, predefinedlist, parentfig)
    % Sample list of items
    items = arrayfun(@(x) sprintf('Item %d',x), 1:20, 'UniformOutput', false);

    % Create UI figure
    fig = uifigure('Name','Item Selection','Position',[100 100 600 400]);

    % Available list
    lblAvail = uilabel(fig,'Text','Available Items','Position',[50 360 120 20]);
    lbAvail = uilistbox(fig,...
        'Position',[50 100 200 260],...
        'Items',items,...
        'Multiselect','on');

    % Selected list
    lblSel = uilabel(fig,'Text','Selected Items','Position',[350 360 120 20]);
    lbSel = uilistbox(fig,...
        'Position',[350 100 200 260],...
        'Multiselect','on');

    % Buttons for moving items
    btnAdd = uibutton(fig,'Text','>',...
        'Position',[270 280 60 30],...
        'ButtonPushedFcn',@(btn,event) moveItems(lbAvail,lbSel));
    btnAddAll = uibutton(fig,'Text','>>',...
        'Position',[270 240 60 30],...
        'ButtonPushedFcn',@(btn,event) moveAll(lbAvail,lbSel));
    btnRemove = uibutton(fig,'Text','<',...
        'Position',[270 200 60 30],...
        'ButtonPushedFcn',@(btn,event) moveItems(lbSel,lbAvail));
    btnRemoveAll = uibutton(fig,'Text','<<',...
        'Position',[270 160 60 30],...
        'ButtonPushedFcn',@(btn,event) moveAll(lbSel,lbAvail));

    % Done button
    uibutton(fig,'Text','Done',...
        'Position',[250 40 100 30],...
        'ButtonPushedFcn',@(btn,event) doneCallback(lbSel));

end

% ---- Helper functions ----
function moveItems(srcList, dstList)
    items = srcList.Items;
    sel   = srcList.Value;

    % Remove from source
    srcList.Items(ismember(items,sel)) = [];

    % Add to destination
    dstList.Items = [dstList.Items(:); sel(:)];
end

function moveAll(srcList, dstList)
    % Move everything from src to dst
    dstList.Items = [dstList.Items(:); srcList.Items(:)];
    srcList.Items = {};
end

function doneCallback(lbSel)
    % Return the selected items
    selectedItems = lbSel.Items;
    disp('Final selected items:');
    disp(selectedItems);
    uialert(lbSel.Parent, sprintf('You selected %d items.', numel(selectedItems)), 'Selection Complete');
end
%}

function mustBeValidItemList(x)
    if ~isempty(x)
        if ~(iscellstr(x) || isstring(x) || iscategorical(x))
            error('Items must be a cell array of char, string array, or 1-D categorical array.');
        end
        if ~isvector(x)
            error('Items must be a 1-D array.');
        end
    end
end