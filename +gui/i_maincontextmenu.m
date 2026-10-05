function cm = i_maincontextmenu(app)
%I_MAINCONTEXTMENU Right-click menu for the main window's cell scatter.
%
%   cm = gui.i_maincontextmenu(app)
%
% Builds the menu once per window and attaches it to the main axes and to
% every scatter drawn in them. The scatter is rebuilt in several places
% (in_UpdateMainPlot, the 2D/3D switch, re-embedding), so rather than
% re-attach it at each, a listener on the axes' ChildAdded event attaches
% it to whatever is added. ChildAdded is undocumented; if it stops firing,
% the menu still opens from the axes background.
%
% Each item runs the menu item it mirrors, through that item's own
% MenuSelectedFcn, so undo snapshots and checks stay in one place. Items
% that act on brushed cells are enabled only while cells are brushed.
% Calling it again is cheap and returns the existing menu.
%
% see also: gui.i_applydisplay

tag = 'scgeatoolMainContextMenu';
cm = findobj(app.UIFigure, 'Type', 'uicontextmenu', 'Tag', tag);
if ~isempty(cm)
    in_attach(app, cm);
    return;
end

cm = uicontextmenu(app.UIFigure, 'Tag', tag);
% Text, the app menu item it runs, and whether it needs brushed cells.
items = {
    'Label Cell Groups...',                    'LabelCellGroupsMenu',                   false
    'Hide Group Labels',                       '',                                      false
    'Show Gene Expression...',                 'SelectedGenesMenu',                     false
    'Expand Brushed Cells to Whole Group',     'ExpandHighlightedCellstoGroupMenu',     true
    'Add Brushed Cells to a New Group...',     'AddBrushedCellstoaNewGroupMenu',        true
    'Merge Groups of Brushed Cells...',        'MergeBrushedCellstoSameGroupMenu',      true
    'Identify Cell Type of Brushed Cells...',  'AnnotateCellTypesforBrushedCellsMenu', true
    'Find Marker Genes for Selected Cells...', 'FindMarkerGenesforBrushedCellsMenu',    true
    'Delete Brushed Cells...',                 'DeleteBrushedCellsMenu',                true
    'Switch Between 2D/3D Embeddings...',      'SwitchBetween2D3DEmbeddingsMenu',       false
    'Refresh Current View',                    'RefreshCurrentViewMenu',                false
    };
separatorsBefore = ["Expand Brushed Cells to Whole Group", "Switch Between 2D/3D Embeddings..."];
handles = gobjects(size(items, 1), 1);
for k = 1:size(items, 1)
    handles(k) = uimenu(cm, 'Text', items{k, 1});
    if any(separatorsBefore == items{k, 1})
        handles(k).Separator = 'on';
    end
    if isempty(items{k, 2})
        handles(k).MenuSelectedFcn = @(~, ~) in_hidelabels(app);
    else
        handles(k).MenuSelectedFcn = @(~, event) in_run(app, items{k, 2}, event);
    end
end
cm.ContextMenuOpeningFcn = @(~, ~) in_refresh(app, handles, items);
cm.UserData = addlistener(app.UIAxes, 'ChildAdded', @(~, event) in_adopt(cm, event.ChildNode));
in_attach(app, cm);
end

function in_attach(app, cm)
app.UIAxes.ContextMenu = cm;
if pkg.i_isvalid(app.h) && isprop(app.h, 'ContextMenu')
    app.h.ContextMenu = cm;
end
end

function in_adopt(cm, child)
% Data tips refuse a custom menu; only chart objects are adopted.
if isvalid(cm) && isprop(child, 'ContextMenu') && ~isa(child, 'matlab.graphics.datatip.DataTip')
    child.ContextMenu = cm;
end
end

function in_run(app, menuName, event)
source = app.(menuName);
if strcmp(source.Enable, 'on')
    source.MenuSelectedFcn(source, event);
end
end

function in_hidelabels(app)
if pkg.i_isvalid(app.h)
    delete(findobj(app.h, 'Type', 'datatip'));
end
end

function in_refresh(app, handles, items)
% Opening is the only moment the brush state and the labels are known.
hasplot = pkg.i_isvalid(app.h);
brushed = hasplot && isprop(app.h, 'BrushData') && any(app.h.BrushData(:));
labelled = hasplot && ~isempty(findobj(app.h, 'Type', 'datatip'));
for k = 1:numel(handles)
    if isempty(items{k, 2})
        ok = labelled;
    else
        ok = strcmp(app.(items{k, 2}).Enable, 'on') && (~items{k, 3} || brushed);
    end
    handles(k).Enable = matlab.lang.OnOffSwitchState(ok);
end
end
