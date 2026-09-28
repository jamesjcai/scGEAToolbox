function callback_CloseAllOthers(src, ~)
% Callback function that closes all other MATLAB figures except the current one.

% Get the handle of the currently active figure
[FigureHandle] = gui.gui_getfigsce(src);

% Visible windows only: findall also returns hidden figures, which are not
% windows the user can see or meant to close.
allFigures = findall(0, 'Type', 'Figure', 'Visible', 'on');
others = allFigures(allFigures ~= FigureHandle);

if isempty(others)
    % Used to do nothing at all, which looked like the menu was broken.
    gui.myHelpdlg(FigureHandle, 'No other figures are open.', '', true);
    return;
end

confirmation = gui.myQuestdlg(FigureHandle, ...
    sprintf('Close %s?', pkg.i_plural(numel(others), 'other figure')), ...
    'Confirmation');
if ~strcmp(confirmation, 'Yes'), return; end

for k = 1:numel(others)
    try
        % CLOSE runs the window's CloseRequestFcn, which may ask first.
        close(others(k));
    catch closeError
        disp(['Failed to close figure: ', closeError.message]);
    end
end

% Counted, not assumed: a window whose own close request was declined,
% or that failed to close, is still open.
left = others(isvalid(others));
if isempty(left)
    gui.myHelpdlg(FigureHandle, 'All other figures have been closed.', '', true);
else
    gui.myHelpdlg(FigureHandle, sprintf('%s still open.', ...
        pkg.i_plural(numel(left), 'figure is', 'figures are')), '', true);
end
end
