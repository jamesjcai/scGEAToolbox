function i_captureframe(fig, tab, filename)
%I_CAPTUREFRAME Write a picture of a figure, or of one tab of it.
%
%   gui.i_captureframe(fig, tab, filename)
%
%   fig       the figure to photograph
%   tab       one of its tabs to photograph instead, [] for the whole figure
%   filename  where to write the picture
%
% EXPORTGRAPHICS draws the graphics themselves; GETFRAME photographs the
% screen. Drawing is better where it works: it needs nothing to be selected or
% even visible, so a window lying over the figure cannot end up in the deck,
% and it renders at print resolution rather than the display's.
%
% A tab is always drawn. That also gets the plots without the tab strip around
% them, and leaves the user's selected tab where they put it.
%
% The picture comes out on a white background whatever the figure looks like
% on screen, because a slide is not a screen. See IN_WHITEBACKDROP: painting
% only the margin white is not enough for a dark figure, so the figure is put
% into the light theme for the duration of the capture and back afterwards. It
% does flash light while that happens; GUI.I_SAVEMAINFIG's PDF and image
% exports have always done the same.
%
% A whole figure is drawn only when EXPORTGRAPHICS can do it without losing
% anything. It cannot draw UI components - a pushbutton, a dropdown, a table -
% and says so with a warning rather than an error, so a figure carrying them
% would otherwise reach the deck with its controls silently missing. Worse
% here than elsewhere: GUI.I_EXPORT2PPTX turns all warnings off around its
% capture loop, so nobody would ever read it. MATLAB's own warning is what
% decides, rather than a list of component classes kept here and guessed at.
%
% see also: gui.i_export2pptx, gui.i_pptxslideplan, gui.i_figtabs

if nargin < 3, filename = ''; end
if nargin < 2, tab = []; end

restorelook = in_whitebackdrop(fig);   %#ok<NASGU>

if ~isempty(tab) && pkg.i_isvalid(tab)
    if in_drawn(tab, filename), return; end

    % EXPORTGRAPHICS refuses a container with nothing drawn in it
    % ('MATLAB:print:EmptyFigureNotSupported'), and there is no reason to
    % assume that is the only thing it will ever refuse. The tab is still on
    % screen, so photograph it there - with the strip around it, which beats
    % a slide that is not in the deck at all.
    in_selecttab(tab);
elseif in_drawn(fig, filename)
    return;
end

frame = getframe(fig);
imwrite(frame.cdata, filename);
end


function tf = in_drawn(h, filename)
% Draw H into FILENAME, and say whether the result can be used. False means
% EXPORTGRAPHICS either refused outright or warned that it left something out,
% and the caller should fall back to photographing the screen.

tf = false;

% LASTWARN is read below, so put back whatever the caller had in it. Anything
% that inspects the last warning after an export is entitled to see its own.
[prevmsg, previd] = lastwarn();
restorelastwarn = onCleanup(@() lastwarn(prevmsg, previd));
lastwarn('', '');

try
    exportgraphics(h, filename, 'Resolution', 150, ...
        'BackgroundColor', 'white');
catch
    % Nothing drawn in it, or any other refusal. The caller has a fallback.
    return;
end

% Still readable under warning('off','all'), which is the whole reason this
% works from inside GUI.I_EXPORT2PPTX's capture loop.
[~, id] = lastwarn();
tf = ~strcmp(id, 'MATLAB:print:ExportappForUIFigureWithUIControl');
end


function c = in_whitebackdrop(fig)
% Put FIG on a white background for the duration of the capture, and give back
% an ONCLEANUP that restores it.
%
% EXPORTGRAPHICS's 'BackgroundColor' paints what is around the plots and
% nothing else. On a dark figure that leaves black axes interiors and pale
% grey titles stranded in a white margin, the titles all but unreadable. What
% decides the axes interiors, the rulers and the text is the theme, so the
% theme is what has to move - the same thing GUI.I_SAVEMAINFIG does before it
% writes a PDF or an image.
%
% The figure's own Color is whitened too. That is for the GETFRAME fallback,
% which photographs the window and would otherwise bring the theme's grey
% along with it; EXPORTGRAPHICS ignores it in favour of the argument above.

undo = {};
if pkg.i_isvalid(fig)
    undo = in_lighten(fig, undo);
end
c = onCleanup(@() in_undoall(undo));
end


function undo = in_lighten(fig, undo)
% Everything worth putting back is read before anything is changed. Setting
% the theme repaints the figure, Color included, so a Color read after that
% line is the light theme's grey and not what the user had - restoring it
% would leave a dark figure looking washed out for the rest of the session.
hascolor = isprop(fig, 'Color');
if hascolor, wascolor = fig.Color; end

if isprop(fig, 'Theme')
    try
        was = fig.Theme.BaseColorStyle;
        if ~strcmp(was, 'light')
            theme(fig, 'light');
            undo{end+1} = @() in_settheme(fig, was);
        end
    catch
        % No THEME on this release, or a figure that will not take one. The
        % background argument still whitens the margin, which is as far as
        % this can go.
    end
end

if hascolor
    fig.Color = [1 1 1];
    % Appended after the theme, and IN_UNDOALL runs them in that order, so a
    % restored theme cannot repaint over the restored colour.
    undo{end+1} = @() in_setcolor(fig, wascolor);
end

% Applied, not merely requested: GETFRAME reads the screen, so the window has
% to have repainted before anything is captured.
drawnow;
end


function in_undoall(undo)
% In the order they were added, not reversed: the theme goes back first and
% the figure's own Color after it, because restoring a theme repaints Color.
for k = 1:numel(undo)
    undo{k}();
end
end


function in_settheme(fig, was)
if ~pkg.i_isvalid(fig), return; end
try
    theme(fig, was);
catch
    % Closed mid-export, or a theme it will no longer accept. Nothing more
    % can be done about the figure's appearance from here.
end
end


function in_setcolor(fig, wascolor)
if pkg.i_isvalid(fig)
    fig.Color = wascolor;
end
end


function in_selecttab(tab)
% Bring TAB forward and let it paint. Without the DRAWNOW, GETFRAME takes the
% screen as it stands and returns whichever tab was already showing.

tg = tab.Parent;
if pkg.i_isvalid(tg) && isprop(tg, 'SelectedTab')
    tg.SelectedTab = tab;
    drawnow;
end
end
