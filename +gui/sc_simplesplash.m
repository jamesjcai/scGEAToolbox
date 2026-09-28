function [fx, v1] = sc_simplesplash(fx, r)
%SC_SIMPLESPLASH Figure-based splash screen, the last fallback of gui.sc_splashscreen.
%   [fx, v1] = gui.sc_simplesplash() shows today's splash picture with the
%   version and a progress bar at 10%, and returns the figure and version.
%   gui.sc_simplesplash(fx, r) moves the progress bar of splash FX to R (0-1).
%   With no pictures in the splash folder it shows nothing and returns [].
%
%   After https://www.mathworks.com/matlabcentral/answers/92259

v1 = '';
if nargin < 2, r = 0.1; end
if nargin < 1 || isempty(fx)
    fx = [];
    splashpng = i_splashimage();
    if splashpng == "", return; end
    im = imread(splashpng);

    fx = figure('MenuBar', 'none', 'ToolBar', 'none', ...
        'Name', '', 'NumberTitle', 'off', 'Color', 'k', ...
        'Position', [0 0 400 300], 'Visible', 'off', ...
        'WindowStyle', 'normal', ...
        'DockControls', 'off', ...
        'Resize', 'off');
    fa = axes('Parent', fx, 'Visible', 'off');
    ih = image(im, 'Parent', fa);
    imxpos = get(ih, 'XData');
    imypos = get(ih, 'YData');
    set(fa, 'Unit', 'Normalized', 'Position', [0, 0, 1, 1]);
    figpos = get(fx, 'Position');
    figpos(3:4) = [imxpos(2) imypos(2)];
    set(fx, 'Position', figpos);

    hold(fa, "on");
    [x, y] = i_progressxy(r);
    plot(fa, x, y, '-', 'LineWidth', 4, 'Color', [0.7 0.7 0.7], 'Tag', 'SplashProgress');
    box(fa, "off");
    axis(fa, "off");
    movegui(fx, 'center');

    text(fa, 20, 50, 'SCGEATOOL', 'Color', 'w', 'FontSize', 16);
    v1 = pkg.i_get_versionnum;
    text(fa, 20, 80, v1, 'Color', [0.7 0.7 0.7], 'FontSize', 12);
    text(fa, 20, 270, 'Loading...', 'Color', [0.7 0.7 0.7], 'FontSize', 12);
    fx.Visible = true;

    % Keep the splash on screen long enough to be seen, even when the app
    % initializes quickly (the caller deletes the figure as soon as startup
    % finishes, which can otherwise make the splash flash by).
    minSplashSeconds = 1.5;
    drawnow;
    pause(minSplashSeconds);
elseif pkg.i_isvalid(fx)
    % Move the one progress line. This used to re-list the splash folder,
    % re-seed the random stream and re-read the picture on every update, and
    % then draw another line on top of the last.
    hLine = findobj(fx, 'Tag', 'SplashProgress');
    if ~isempty(hLine)
        [x, y] = i_progressxy(r);
        set(hLine(1), 'XData', x, 'YData', y);
    end
end
end

function [x, y] = i_progressxy(r)
% Progress-bar line for ratio R along the bottom of the 400-pixel-wide picture.
X = 10:390;
x = X(1:round(length(X)*min(max(r, 0), 1)));
y = 290*ones(size(x));
end
