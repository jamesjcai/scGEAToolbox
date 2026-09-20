function [out] = e_makedisplaystruct

% Look-and-feel of the main scatter, carried on the SCE so it survives a
% save and reload. Empty means "never set", and the app then draws with
% its own defaults, so an SCE made before this existed behaves as before.
%
% Marker and SizeData travel together: a point marker needs a far larger
% size than a shaped one to stay visible, so restoring one without the
% other gives an invisible or blotted plot.
%
% Colormap is the N-by-3 matrix rather than a name because
% GUI.CALLBACK_PICKCOLORMAP applies matrices, several of which have no
% name to record.
%
% View is [azimuth elevation]. It is also what tells 2D from 3D here: the
% app calls a plot 2D when the elevation is 90, looking straight down, so
% storing the angle stores the 2D/3D choice with it.

out = struct('Marker', [], 'SizeData', [], 'Colormap', [], ...
    'View', [], 'Version', 1);
out = orderfields(out);
