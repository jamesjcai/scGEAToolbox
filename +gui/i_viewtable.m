function hFig = i_viewtable(T, ParentFigureHandle, tabTitles)
%I_VIEWTABLE Show one or more tables in a window.
%   hFig = gui.i_viewtable(T, parentfig) shows table T.
%   hFig = gui.i_viewtable({T1, T2}, parentfig, ["Title 1", "Title 2"])
%   shows each table on its own tab of one window, rather than a window
%   per table stacked over the main app.

if nargin < 3, tabTitles = []; end
if nargin < 2, ParentFigureHandle = []; end

% https://www.mathworks.com/matlabcentral/answers/254690-how-can-i-display-a-matlab-table-in-a-figure

hFig = uifigure('Visible', false);
if istable(T)
    in_addtable(hFig, T);
else
    tabgp = uitabgroup(hFig, 'Units', 'normalized', 'Position', [0, 0, 1, 1]);
    for k = 1:numel(T)
        in_addtable(uitab(tabgp, 'Title', tabTitles(k)), T{k});
    end
end
gui.i_movegui2parent(hFig, ParentFigureHandle);
set(hFig, 'visible', 'on');
end


function in_addtable(parent, T)
uitable(parent, 'Data', T, 'ColumnName', T.Properties.VariableNames, ...
'RowName', T.Properties.RowNames, 'Units', ...
'Normalized', 'Position', [0.05, 0.05, 0.92, 0.90]);
end
