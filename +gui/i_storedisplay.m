function i_storedisplay(app)
%I_STOREDISPLAY Record the scatter's current look on the SCE.
%
%   gui.i_storedisplay(app)
%
%   Called whenever the user changes the marker or the colormap, so the
%   SCE always carries the look currently on screen. Writing here rather
%   than at save time is what makes every save path keep it - the Save
%   dialog, Export/Save Data, and a plain "save f sce" at the prompt all
%   just write the object.
%
%   Anything that is not an app is ignored: the pick callbacks are also
%   reachable with a bare figure, where there is no SCE to write to.
%
% See also GUI.I_APPLYDISPLAY, PKG.E_MAKEDISPLAYSTRUCT.

if nargin < 1 || isempty(app) || ~isa(app, 'matlab.apps.AppBase')
    return;
end
if ~isprop(app, 'sce') || isempty(app.sce) || app.sce.NumCells == 0
    return;
end

d = app.sce.struct_display;
if ~isstruct(d) || ~all(isfield(d, {'Marker', 'SizeData', 'Colormap'}))
    d = pkg.e_makedisplaystruct;
end

if isprop(app, 'h') && ~isempty(app.h) && all(pkg.i_isvalid(app.h))
    h1 = app.h(1);
    d.Marker = h1.Marker;
    if isprop(h1, 'SizeData')
        d.SizeData = h1.SizeData;
    end
end

try
    d.Colormap = colormap(app.UIAxes);
    [az, el] = view(app.UIAxes);
    d.View = [az, el];
catch
    % No axes yet, or none with a colormap; leave whatever was stored.
end

app.sce.struct_display = d;
end
