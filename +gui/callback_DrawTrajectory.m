function f_traj = callback_DrawTrajectory(src, switchTo2DFn, refreshFn)
% callback_DrawTrajectory  Draw/compute a pseudotime trajectory curve.
%
%   f_traj = callback_DrawTrajectory(src)
%   f_traj = callback_DrawTrajectory(src, switchTo2DFn, refreshFn)
%
%   src          - App handle (matlab.apps.AppBase) or legacy figure handle
%   switchTo2DFn - (optional) function handle: @() PushTool_Switch2D3DEmbeddingsClicked(app,[])
%                  Called when the user asks to switch from 3D to 2D before manual drawing.
%   refreshFn    - (optional) function handle: @(keepview,keepcolr) in_RefreshAll(app,kv,kc)
%                  Called when removing an existing trajectory before redrawing.
%
%   f_traj       - Graphics handle of the plotted trajectory line, or [] if cancelled.

f_traj = [];

if nargin < 2, switchTo2DFn = []; end
if nargin < 3, refreshFn    = []; end

[FigureHandle, sce] = gui.gui_getfigsce(src);

% Extract app-specific handles when src is an App Designer app
if isa(src, 'matlab.apps.AppBase')
    UIAxes     = src.UIAxes;
    cur_f_traj = src.f_traj;
    scatter_h  = src.h;
else
    UIAxes     = [];
    cur_f_traj = [];
    scatter_h  = [];
end

% --- Preamble warning ---
if ~strcmpi('Yes', gui.myQuestdlg(FigureHandle, ['This function is ' ...
        'not intended for use with t-SNE or UMAP embeddings, since ' ...
        'those methods tend to break up the data into separate ' ...
        'clusters, which can hide the smooth transitions you''d ' ...
        'expect in developmental trajectories. Use PHATE ' ...
        'embedding instead. Continue?']))
    return;
end

% --- Remove existing trajectory if present ---
if ~isempty(cur_f_traj) && isvalid(cur_f_traj) && isgraphics(cur_f_traj, 'line')
    if strcmp({cur_f_traj.Visible}, 'on')
        switch gui.myQuestdlg(FigureHandle, 'Remove existing trajectory curve?', '')
            case 'Yes'
                if ~isempty(refreshFn)
                    refreshFn(true, true);  % in_RefreshAll(app, keepview=true, keepcolr=true)
                end
            case 'No'
            otherwise
                return;
        end
    end
end

% --- Choose method ---
if license('test', 'curve_fitting_toolbox') && ~isempty(which('cscvn'))
    answer = gui.myQuestdlg(FigureHandle, 'Which method?', '', ...
        {'splinefit', 'princurve', 'manual'}, 'splinefit');
else
    answer = gui.myQuestdlg(FigureHandle, 'Which method?', '', ...
        {'splinefit', 'princurve'}, 'splinefit');
end

justload = false;

switch answer
    case 'splinefit'
        dim = 1;
        [t, xyz1] = pkg.i_pseudotime_by_splinefit(sce.s, dim, false);
        pseudotimemethod = 'splinefit';

    case 'princurve'
        try
            [t, xyz1] = pkg.i_pseudotime_by_princurve(sce.s, false);
            pseudotimemethod = 'princurve';
        catch ME
            gui.myErrordlg(FigureHandle, "Runtime error.");
            return;
        end

    case 'manual'
        if license('test', 'curve_fitting_toolbox') && ~isempty(which('cscvn'))
            % If current embedding is 3D, offer to switch to 2D first
            if ~isempty(scatter_h) && isvalid(scatter_h) && ~isempty(scatter_h.ZData)
                switch gui.myQuestdlg(FigureHandle, ...
                        ['This function does not work for 3D ' ...
                        'embedding. Continue to switch to 2D?'])
                    case 'Yes'
                        if ~isempty(switchTo2DFn)
                            switchTo2DFn();
                        end
                    otherwise
                        return;
                end
            end
            % Re-check: if still 3D after switch attempt, bail
            if ~isempty(scatter_h) && isvalid(scatter_h) && ~isempty(scatter_h.ZData), return; end

            answer2 = gui.myQuestdlg(FigureHandle, ...
                'Draw trajectory curve or load saved curve and pseudotime?', ...
                '', {'Draw Curve', 'Load Saved', 'Cancel'}, 'Draw Curve');

            switch answer2
                case 'Load Saved'
                    [file, path] = uigetfile('*.mat', 'Select a MAT-file to Load');
                    if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure')
                        figure(FigureHandle);
                    end
                    if isequal(file, 0)
                        disp('User canceled the file selection.');
                        return;
                    end
                    fullFileName = fullfile(path, file);
                    loadedData = load(fullFileName);
                    if isfield(loadedData, 't') && isfield(loadedData, 'xyz1') ...
                            && isfield(loadedData, 'pseudotimemethod')
                        t               = loadedData.t;
                        xyz1            = loadedData.xyz1;
                        pseudotimemethod = loadedData.pseudotimemethod;
                    else
                        gui.myErrordlg(FigureHandle, 'Not a valid .mat file.', '');
                        return;
                    end
                    justload = true;

                case 'Draw Curve'
                    if isempty(UIAxes)
                        gui.myErrordlg(FigureHandle, 'Manual drawing requires an App Designer figure.', '');
                        return;
                    end
                    fx   = gui.myFigure(FigureHandle);
                    fig2 = fx.FigHandle;
                    ax2  = fx.AxHandle;
                    gui.copyUIAxesToFigure(UIAxes, ax2);
                    fx.show(FigureHandle);
                    fig2.WindowStyle = "modal";
                    hold(ax2, "on");
                    x = []; y = [];
                    while true
                        [xi, yi, button] = ginput(1);
                        if isempty(button) || button == 13
                            break;
                        end
                        x = [x; xi];
                        y = [y; yi];
                        plot(ax2, xi, yi, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                    end
                    pause(1);
                    fx.closeFigure;

                    hold(UIAxes, "on");
                    for k = 1:length(x)
                        plot(UIAxes, x(k), y(k), 'ro', 'MarkerSize', 8, 'LineWidth', 2);
                    end

                    splineCurve = cscvn([x'; y']);
                    [xyz1] = fnplt(splineCurve);
                    xyz1 = xyz1';

                    [t] = dsearchn(xyz1, sce.s);
                    t   = (t + randn(size(t)))';
                    t   = normalize(t, 'range');
                    t   = t(:);
                    pseudotimemethod = 'manual';

                otherwise
                    return;
            end
        else
            return;
        end

    otherwise
        return;
end

% --- Plot trajectory ---
if isempty(UIAxes), return; end  % can't draw without axes

% Validate dimensions
if size(xyz1, 2) < 2, return; end

% Interpolate if needed (3D case only)
if size(xyz1, 2) >= 3 && size(xyz1, 1) < 0.8 * sce.NumCells
    xyz1 = pkg.i_interp3d(xyz1, sce.NumCells);
end

hold(UIAxes, 'on');

if size(xyz1, 2) >= 3
    f_traj    = plot3(UIAxes, xyz1(:, 1), xyz1(:, 2), xyz1(:, 3), '-r', 'linewidth', 2);
    start_pos = [xyz1(1, 1), xyz1(1, 2), xyz1(1, 3)];
    end_pos   = [xyz1(end, 1), xyz1(end, 2), xyz1(end, 3)];
else
    f_traj    = plot(UIAxes, xyz1(:, 1), xyz1(:, 2), '-r', 'linewidth', 2);
    start_pos = [xyz1(1, 1), xyz1(1, 2)];
    end_pos   = [xyz1(end, 1), xyz1(end, 2)];
end

% Determine text properties based on MATLAB version and theme
if isMATLABReleaseOlderThan('R2025a')
    text_props = {'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k'};
else
    style = FigureHandle.Theme.BaseColorStyle;
    if strcmp(style, "light")
        text_props = {'fontsize', 10, 'FontWeight', 'bold', 'BackgroundColor', 'w', 'EdgeColor', 'k'};
    else % dark theme
        text_props = {'fontsize', 10, 'FontWeight', 'bold', 'EdgeColor', 'w'};
    end
end

t1 = text(UIAxes, start_pos(1), start_pos(2), start_pos(end), 'Start', text_props{:});
t2 = text(UIAxes, end_pos(1),   end_pos(2),   end_pos(end),   'End',   text_props{:});

hold(UIAxes, 'off');
drawnow;
pause(2);

% --- Swap start/end? ---
if ~strcmp(answer, 'manual')
    switch gui.myQuestdlg(FigureHandle, 'Swap ''Start'' and ''End''?', '')
        case 'Yes'
            t1.String = 'End';
            t2.String = 'Start';
            t = 1 - t;
        case 'Cancel'
            return;
    end
end

% --- Save pseudotime to sce ---
tag = sprintf('%s_pseudotime', pseudotimemethod);
try
    idx = find(contains(sce.list_cell_attributes(1:2:end), tag));
catch ME
    idx = [];
    warning(ME.message);
end
if ~isempty(idx)
    sce.list_cell_attributes{idx*2} = t;
    fprintf('%s is updated.\n', upper(tag));
else
    sce.list_cell_attributes{end+1} = tag;
    sce.list_cell_attributes{end+1} = t;
    fprintf('%s is saved.\n', upper(tag));
end

% --- Optionally save manual curve to .mat ---
if strcmp(answer, 'manual') && ~justload
    switch gui.myQuestdlg(FigureHandle, ...
            'Save manual trajectory curve and pseudotime to an .mat file?', '')
        case 'Yes'
            [file, path] = uiputfile('*.mat', 'Save as', 'pseudotime_manual_trajectory.mat');
            if isvalid(FigureHandle) && isa(FigureHandle, 'matlab.ui.Figure')
                figure(FigureHandle);
            end
            if isequal(file, 0) || isequal(path, 0)
                disp('User canceled the file selection.');
                return;
            end
            fullFileName = fullfile(path, file);
            save(fullFileName, 'xyz1', 't', 'pseudotimemethod');
            disp(['Variables saved to ', fullFileName]);
        case 'No'
        case 'Cancel'
            return;
        otherwise
            return;
    end
end

% --- View expression of selected genes ---
switch gui.myQuestdlg(FigureHandle, 'View expression of selected genes', '')
    case 'Yes'
        gui.sc_pseudotimegenes(sce, t, FigureHandle);
end
end
