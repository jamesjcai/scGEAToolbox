function i_viewtable(T)

%T.Properties.RowNames = T.hvg;
%T = T(:, 2:end);

% https://www.mathworks.com/matlabcentral/answers/254690-how-can-i-display-a-matlab-table-in-a-figure

% LastName = {'Smith';'Johnson';'Williams';'Jones';'Brown'};
% Age = [38;43;38;40;49];
% Height = [71;69;64;67;64];
% Weight = [176;163;131;133;119];
% T = table(Age,Height,Weight,'RowNames',LastName);

hFigure = uifigure('Visible', false);
%set(hFigure, 'MenuBar', 'none');
%set(hFigure, 'ToolBar', 'none');
uitable(hFigure, 'Data', T{:, :}, 'ColumnName', T.Properties.VariableNames, ...
    'RowName', T.Properties.RowNames, 'Units', ...
    'Normalized', 'Position', [0.03, 0.03, 0.94, 0.93]);
movegui(hFigure, 'center');
set(hFigure, 'visible', 'on');
