function [pt] = i_3dcamera(tb, prefix, flatview, parentfig, parentax)
%I_3DCAMERA Toolbar button that records a rotating video of a 3-D plot.
%   pt = gui.i_3dcamera(tb, prefix, flatview, parentfig, parentax) adds a
%   push tool to toolbar TB. PREFIX starts the video file name, FLATVIEW
%   picks the camera path (true: a lower, flatter orbit), PARENTFIG is the
%   window dialogs are centred on and PARENTAX the axes rotated -- when it is
%   empty, the current axes at the moment of the click.
%
%   This replaces gui.gui_3dcamera, a near-copy without PARENTAX. This one
%   used to ignore FLATVIEW, overwriting it with RAND>0.5 on every click (so
%   the camera path was a coin flip, drawn from the global random stream),
%   and it opened the result with WINOPEN, which exists only on Windows.

if nargin < 5, parentax = []; end
if nargin < 4, parentfig = []; end
if nargin < 3 || isempty(flatview), flatview = false; end
if nargin < 2, prefix = ''; end
if nargin < 1
    hFig = gcf;
    tb = uitoolbar('Parent', hFig);
end
% Raising is about keeping focus on a window the user is already looking
% at. FIGURE() also forces Visible on, so doing it to a hidden figure
% shows it: GUI.MYFIGURE builds its toolbar before the plot is drawn, and
% this popped that half-built figure up as an empty window that vanished
% again when SHOW finally positioned it. Only raise what is already up.
% == "on" rather than strcmp: Visible is a matlab.lang.OnOffSwitchState,
% which strcmp never matches against a char.
if ~isempty(parentfig) && pkg.i_isvalid(parentfig) && parentfig.Visible == "on"
    figure(parentfig);
    cleanupObj = onCleanup(@() gui.i_raisefig(parentfig));
end
pt = uipushtool(tb, 'Separator', 'on');

try
    mfolder = fileparts(mfilename('fullpath'));
    [ptImage, map] = imread(fullfile(mfolder, '..', 'assets', 'Images', 'camera.jpg'));
    if ~isempty(map), ptImage = ind2rgb(ptImage, map); end
catch
    ptImage = rand(16, 16, 3);
end
pt.CData = ptImage;
pt.Tooltip = 'Make video snapshot';
pt.ClickedCallback = @camera3dmp4;


    function camera3dmp4(~, ~)
        answer = gui.myQuestdlg(parentfig, 'Make video snapshot?');
        if ~strcmp(answer, 'Yes'), return; end

        ax = parentax;
        if isempty(ax) || ~isvalid(ax)
            ax = gca;
        end
        fig = ancestor(ax, 'figure');

        [caz, cel] = view(ax);
        OptionZ.FrameRate = 15;
        OptionZ.Duration = 5.5;
        OptionZ.Periodic = true;
        fname = tempname;
        if ~isempty(prefix)
            [a1, b1] = fileparts(fname);
            b1 = sprintf('%s_%s', prefix, b1);
            fname = fullfile(a1, b1);
        end

        % Capture and restore rather than a bare off/on pair: the pair leaves
        % warnings disabled for the rest of the session if anything between the
        % two lines throws, and its 'on' re-enables warnings the caller may have
        % silenced deliberately instead of restoring what they had.
        warnState = warning();
        restoreWarn = onCleanup(@() warning(warnState));
        warning('off', 'all');
        if flatview
            views = [-20, 50; -110, 65; -190, 80; -290, 60; -380, 40];
        else
            views = [-20, 10; -110, 10; -190, 80; -290, 10; -380, 10];
        end
        try
            i_captureFigVid(views, fname, OptionZ, fig, ax);
        catch
            % video export is optional; restore view and continue if codec missing
        end
        view(ax, caz, cel);

        vfile = '';
        for ext = [".mp4", ".avi"]
            if isfile(fname + ext)
                vfile = char(fname + ext);
                break;
            end
        end
        if isempty(vfile)
            gui.myWarndlg(parentfig, 'No video was written (video codec not available?).');
        elseif ispc
            winopen(tempdir);
            pause(1);
            winopen(vfile);
        else
            gui.myHelpdlg(parentfig, sprintf('Video saved to %s', vfile));
        end
    end

end


function i_captureFigVid(ViewZ, FileName, OptionZ, fig, ax)
% Record the axes AX rotating through the view angles VIEWZ (rows of
% [azimuth, elevation]) to FILENAME. After CaptureFigVid by Alan Jennings
% (Air Force Institute of Technology): OptionZ.FrameRate, .Duration (spaces
% the views over that many seconds) and .Periodic (drops the final view so
% the video loops cleanly). MPEG-4 is written on Windows, the VideoWriter
% default elsewhere.
if nargin < 3
    OptionZ = struct([]);
end

% check orientation of ViewZ, should be two columns and >=2 rows
if size(ViewZ, 2) > size(ViewZ, 1)
    ViewZ = ViewZ.';
end
if size(ViewZ, 2) > 2
    warning('AJennings:VidWrite', ...
        'Views should have n rows and only 2 columns. Deleting extraneous input.');
    ViewZ = ViewZ(:, 1:2);
end

if ispc
    daObj = VideoWriter(FileName, 'MPEG-4');
else
    daObj = VideoWriter(FileName);
end
if isfield(OptionZ, 'FrameRate')
    daObj.FrameRate = OptionZ.FrameRate;
end
if isfield(OptionZ, 'Duration') % space out view angles
    temp_n = round(OptionZ.Duration*daObj.FrameRate); % number frames
    temp_p = (temp_n - 1) / (size(ViewZ, 1) - 1); % length of each interval
    ViewZ_new = zeros(temp_n, 2);
    for inis = 1:(size(ViewZ, 1) - 1)
        ViewZ_new(round(temp_p*(inis - 1)+1):round(temp_p*inis+1), :) = ...
            [linspace(ViewZ(inis, 1), ViewZ(inis+1, 1), ...
            round(temp_p*inis)-round(temp_p*(inis - 1))+1).', ...
            linspace(ViewZ(inis, 2), ViewZ(inis+1, 2), ...
            round(temp_p*inis)-round(temp_p*(inis - 1))+1).'];
    end
    ViewZ = ViewZ_new;
end
if length(ViewZ) == 2 % only initial and final given
    ViewZ = [linspace(ViewZ(1, 1), ViewZ(end, 1)).', ...
        linspace(ViewZ(1, 2), ViewZ(end, 2)).'];
end
if isfield(OptionZ, 'Periodic') && OptionZ.Periodic
    ViewZ = ViewZ(1:(end -1), :); % remove last sample
end
open(daObj);
for kathy = 1:size(ViewZ, 1)
    view(ax, ViewZ(kathy, :));
    drawnow;
    writeVideo(daObj, getframe(fig)); % the figure: the axes change size with the view
end
close(daObj);
end
