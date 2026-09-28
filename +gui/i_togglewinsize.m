function i_togglewinsize(src, ~, hFig, frac)
%I_TOGGLEWINSIZE Toolbar callback: switch a plot window between two sizes.
%
%   Every gui.myFigure window already has this button in its standard
%   toolbar (hx.SizeButton). Other figures add it themselves:
%
%   gui.i_addbutton2fig(tb, 'off', {@gui.i_togglewinsize, hFig}, ...
%       'scale-frame-enlarge.jpg', 'Enlarge Window');
%
%   The first click centres the window at FRAC (default 0.8) of the monitor
%   it is on; the next click puts it back exactly where and how big it was.
%   Maximizing was tried first and filled the screen, more than a plot needs.
%
%   The position to restore is kept in the button's UserData, so the helper
%   holds no state of its own and any number of windows can use it at once.
%   The button's icon and tooltip follow the state, so a click always does
%   what the tooltip says.

if nargin < 4 || isempty(frac), frac = 0.8; end
if isempty(hFig) || ~isvalid(hFig), return; end
if strcmp(hFig.WindowStyle, 'docked')
    % A docked window's size belongs to the desktop, not to the figure.
    gui.myHelpdlg(hFig, 'Undock the figure to resize it.', '');
    return;
end

hFig.Units = 'pixels';
if isempty(src.UserData)
    smallPos = hFig.Position;
    monitors = get(groot, 'MonitorPositions');
    cx = smallPos(1) + smallPos(3)/2;
    cy = smallPos(2) + smallPos(4)/2;
    inMon = cx >= monitors(:, 1) & cx < monitors(:, 1) + monitors(:, 3) & ...
        cy >= monitors(:, 2) & cy < monitors(:, 2) + monitors(:, 4);
    m = find(inMon, 1);
    if isempty(m), m = 1; end   % window off-screen: use the primary
    mon = monitors(m, :);
    wh = round(frac*mon(3:4));
    hFig.Position = [mon(1:2) + round((mon(3:4) - wh)/2), wh];
    src.UserData = smallPos;
    in_seticon(src, 'scale-frame-reduce.jpg', 'Restore Window Size');
else
    hFig.Position = src.UserData;
    src.UserData = [];
    in_seticon(src, 'scale-frame-enlarge.jpg', 'Enlarge Window');
end
end

function in_seticon(btn, imgfile, tip)
imgPath = fullfile(fileparts(mfilename('fullpath')), '..', ...
    'assets', 'Images', imgfile);
try
    [img, map] = imread(imgPath);
    if ~isempty(map), img = ind2rgb(img, map); end
    btn.CData = img;
catch
    % Icon unavailable; the tooltip still says what a click does.
end
btn.Tooltip = tip;
end
