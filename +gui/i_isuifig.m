function y = i_isuifig(FigureHandle)
% https://www.mathworks.com/matlabcentral/answers/348387-distinguish-uifigure-from-figure-programmatically?utm_source=chatgpt.com
y = false;
if isempty(FigureHandle), return; end

y = isprop(FigureHandle,'isUIFigure');

% y = matlab.ui.internal.isUIFigure(FigureHandle);

% if ~isempty(f) && isempty(get(f,'JavaFrame_I'))
%     bool = true;
% else
%     bool = false;
% end

end

%{
fig = uifigure('Name', 'Sample App', 'Position', [100 100 400 300]);
selection = uiconfirm(fig, 'Do you want to continue?', 'Confirm', ...
                      'Options', {'Yes', 'No'}, ...
                      'DefaultOption', 1, ...
                      'CancelOption', 2);

fig = uifigure('Name', 'Sample App', 'Position', [100 100 400 300]);
d = uiprogressdlg(fig, 'Title', 'Please Wait', 'Message', 'Loading...');
pause(2); % Simulate a task
d.Value = 0.5;
pause(2); % Simulate a task
d.Value = 1.0;
close(d);

fig = uifigure('Name', 'Sample App', 'Position', [100 100 400 300]);
uialert(fig, 'This is an informational message.', 'Information', 'Icon', 'info');


fig = uifigure('Name', 'Sample App', 'Position', [100 100 400 300]);
selection = uiconfirm(fig, 'Do you want to continue?', 'Confirm', ...
                      'Options', {'Yes', 'No', '33'}, ...
                      'DefaultOption', 1, ...
                      'CancelOption', 2)


%}