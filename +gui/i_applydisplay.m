function i_applydisplay(app, applyview)
%I_APPLYDISPLAY Put the SCE's stored look back on the scatter.
%
%   gui.i_applydisplay(app)
%   gui.i_applydisplay(app, applyview)
%
%   Called after the scatter is drawn or redrawn. An SCE with nothing
%   stored - one made before the setting existed, or never customised -
%   leaves every default alone, so this is a no-op on anything the user
%   has not deliberately changed.
%
%   APPLYVIEW (default true) governs the camera angle alone. A redraw that
%   is preserving what is on screen passes false: the stored angle is only
%   refreshed when something is picked or the file is saved, so a plot the
%   user had just rotated by hand would otherwise snap back to the last
%   stored angle on the next redraw.
%
%   It is also where the scatter's data tips are labelled, since this is
%   the one call every redraw path makes after the points are drawn.
%
% See also GUI.I_STOREDISPLAY, GUI.I_CELLTYPEDATATIP, PKG.E_MAKEDISPLAYSTRUCT.

if nargin < 2 || isempty(applyview), applyview = true; end
if nargin < 1 || isempty(app) || ~isa(app, 'matlab.apps.AppBase')
    return;
end
if ~isprop(app, 'sce') || isempty(app.sce)
    return;
end

hasscatter = isprop(app, 'h') && ~isempty(app.h) && all(pkg.i_isvalid(app.h));

% Before the struct_display guard below: a data tip that names the cell is
% wanted on every SCE, not only on one that has a stored look.
if hasscatter
    gui.i_celltypedatatip(app.h, app.sce);
end

d = app.sce.struct_display;
if ~isstruct(d)
    return;
end

if hasscatter && isfield(d, 'Marker') && ~isempty(d.Marker)
    set(app.h, 'Marker', d.Marker);
end
% Size goes with the marker: a point needs a much bigger SizeData than a
% shape, so restoring one without the other leaves the plot unreadable.
if hasscatter && isfield(d, 'SizeData') && ~isempty(d.SizeData) ...
        && isprop(app.h(1), 'SizeData')
    set(app.h, 'SizeData', d.SizeData);
end

% A grouped scatter indexes the colormap one row per group - both
% GUI.I_GSCATTER3 and GUI.CALLBACK_PICKCOLORMAP build it that size - so a
% colormap stored under a different grouping no longer fits the plot.
% Restoring it anyway is what flattened a twelve-group plot to one colour:
% the SCE had been looked at ungrouped, storing a single-row map, and
% every later redraw put that one row back over the twelve the redraw had
% just set.
if hasscatter && isfield(d, 'Colormap') && ~isempty(d.Colormap) ...
        && size(d.Colormap, 2) == 3 && i_colormapfits(app.h(1), d.Colormap)
    try
        colormap(app.UIAxes, d.Colormap);
    catch
        % A stored colormap that this release rejects is not worth
        % failing a redraw over; the default one stays.
    end
end

% Only a plot with real depth gets its angle back. Restoring a tilted view
% onto a flat embedding would tip a 2D scatter on its side, and a flat one
% has only the one view worth having anyway.
if applyview && hasscatter && isfield(d, 'View') && numel(d.View) == 2
    z = app.h(1).ZData;
    if ~isempty(z) && ~isscalar(unique(z))
        try
            view(app.UIAxes, d.View(1), d.View(2));
        catch
            % Leave whatever the caller set.
        end
    end
end
end

function tf = i_colormapfits(h, cmap)
%I_COLORMAPFITS Whether CMAP can still colour what H is plotting.
%
%   A grouped scatter carries CData 1..k, one group per row of the
%   colormap, so a stored map of a different height merges groups into one
%   colour or spreads them across the wrong rows. Anything else - a
%   continuous score, a truecolor CData, a flat colour - is drawn by
%   interpolating the map across the axes' CLim, where any height works
%   and the stored choice is still the user's.
tf = true;
if isempty(h) || ~all(pkg.i_isvalid(h)) || ~isprop(h, 'CData'), return; end

cdata = h.CData;
if isempty(cdata) || ~isvector(cdata) || size(cdata, 2) == 3, return; end
cdata = full(double(cdata(:)));
if any(~isfinite(cdata)), return; end

% Grouped means FINDGROUPS output: every integer from 1 to k, all present.
u = unique(cdata);
if ~isequal(u(:)', 1:numel(u)), return; end

tf = size(cmap, 1) == numel(u);
end
