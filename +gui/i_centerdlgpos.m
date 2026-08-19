function pos = i_centerdlgpos(parentfig, dlgSize)
% Returns [x, y, w, h] centered on parentfig and clamped to its monitor.
%
% Monitor detection strategy:
%   1. Area-overlap on parentfig.Position (works for non-maximized windows).
%   2. Mouse-cursor fallback when the parent is maximized: on secondary
%      monitors with DPI different from the primary, MATLAB reports the
%      pre-maximize restore position rather than the actual on-screen bounds,
%      so Position-based detection picks the wrong monitor. The cursor is
%      almost always on the same monitor as the app when a dialog is opened.
w = dlgSize(1);
h = dlgSize(2);
monitors = get(groot, 'MonitorPositions');

if isempty(parentfig) || ~pkg.i_isvalid(parentfig)
    ss = get(groot, 'ScreenSize');
    pos = [ss(1)+(ss(3)-w)/2, ss(2)+(ss(4)-h)/2, w, h];
    return;
end

pPos = parentfig.Position;

% Area-overlap: pick the monitor that contains the most of parentfig.
xOvlp = max(0, min(pPos(1)+pPos(3), monitors(:,1)+monitors(:,3)) - max(pPos(1), monitors(:,1)));
yOvlp = max(0, min(pPos(2)+pPos(4), monitors(:,2)+monitors(:,4)) - max(pPos(2), monitors(:,2)));
[maxArea, idx] = max(xOvlp .* yOvlp);

% Fall back to mouse cursor when Position is unreliable:
%   - maxArea == 0: parentfig is entirely off every known monitor.
%   - WindowState == 'maximized': Position may show restore bounds, not the
%     actual monitor (MATLAB limitation on cross-DPI multi-monitor setups).
useMouseFallback = (maxArea == 0);
if ~useMouseFallback
    try
        useMouseFallback = strcmp(parentfig.WindowState, 'maximized');
    catch
        % older releases or non-figure parents have no WindowState; leave default
    end
end
if useMouseFallback
    idx = i_cursor_monitor_idx(monitors, idx);
end

mon = monitors(idx, :);

% Center dialog on the parent; clamp with a small inset so the dialog
% never lands on a monitor seam (MATLAB/Windows can snap or hide a window
% placed exactly at the pixel boundary between two monitors).
margin = 10;
x = pPos(1) + (pPos(3)-w)/2;
y = pPos(2) + (pPos(4)-h)/2;
x = max(mon(1)+margin, min(x, mon(1)+mon(3)-w-margin));
y = max(mon(2)+margin, min(y, mon(2)+mon(4)-h-margin));

% fprintf('[i_centerdlgpos] pPos=[%g %g %g %g] WindowState=%s useMouse=%d mon=[%g %g %g %g] pos=[%g %g %g %g]\n', ...
%     pPos(1), pPos(2), pPos(3), pPos(4), ...
%     parentfig.WindowState, useMouseFallback, ...
%     mon(1), mon(2), mon(3), mon(4), x, y, w, h);

pos = round([x, y, w, h]);
end

% -------------------------------------------------------------------------
function idx = i_cursor_monitor_idx(monitors, defaultIdx)
% Return the index into monitors of the monitor under the mouse cursor.
% groot's PointerLocation is already reported in the same coordinate space
% as MonitorPositions, so no Java AWT / DPI bridging is needed.
idx = defaultIdx;
try
    pl = get(groot, 'PointerLocation'); % [x y] in MATLAB screen pixels
    mx = pl(1);
    my = pl(2);

    % Find which MATLAB monitor contains the mouse.
    inMon = mx >= monitors(:,1) & mx < monitors(:,1)+monitors(:,3) & ...
            my >= monitors(:,2) & my < monitors(:,2)+monitors(:,4);
    k = find(inMon, 1);
    if ~isempty(k)
        idx = k;
    end
catch ME
    fprintf('[i_centerdlgpos] cursor fallback error: %s\n', ME.message);
end
end
