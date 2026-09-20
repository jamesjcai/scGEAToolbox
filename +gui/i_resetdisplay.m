function i_resetdisplay(app)
%I_RESETDISPLAY Forget the look stored on the SCE, back to app defaults.
%
%   gui.i_resetdisplay(app)
%
%   Called by the two Refresh entry points before they redraw. Clearing the
%   struct is the whole reset: with every field empty GUI.I_APPLYDISPLAY is
%   a no-op, so the marker, size, and colormap that GUI.I_GSCATTER3 sets on
%   the redraw are what stays on screen. An SCE that never had a look picked
%   is already in this state, so a Refresh on one changes nothing.
%
%   This is deliberately not a plain Refresh: every other caller of
%   IN_REFRESHALL is finishing some other operation - a reclustering, an
%   annotation, an undo - and has no business discarding a colormap the user
%   picked. Only the Refresh button and the Refresh Current View menu item
%   reset.
%
%   A snapshot is taken first, so a Refresh aimed at the wrong window is
%   recoverable from Edit > Undo like any other operation. It has to be
%   taken here rather than in the two handlers: STRUCT_DISPLAY is mutated on
%   a handle object, so by the time the assignment below has run there is
%   nothing left to copy - see GUI.I_SNAPSHOT.
%
%   A Refresh with nothing to discard returns before that. Undo holds one
%   level, and overwriting it for a Refresh that changed nothing would cost
%   the user the cell deletion or merge sitting behind it.
%
% See also GUI.I_STOREDISPLAY, GUI.I_APPLYDISPLAY, GUI.I_SNAPSHOT,
% PKG.E_MAKEDISPLAYSTRUCT.

if nargin < 1 || isempty(app) || ~isa(app, 'matlab.apps.AppBase')
    return;
end
if ~isprop(app, 'sce') || isempty(app.sce)
    return;
end

cleared = pkg.e_makedisplaystruct;
if isequal(app.sce.struct_display, cleared)
    return;
end

gui.i_snapshot(app, 'Display Reset');
app.sce.struct_display = cleared;
end
