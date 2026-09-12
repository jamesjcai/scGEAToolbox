function existing = i_dlgregister(parentfig, d)
%I_DLGREGISTER  At most one picker dialog per parent figure.
%
%   existing = GUI.I_DLGREGISTER(parentfig) returns the picker dialog
%   already open on PARENTFIG, raising it to the front, or [] if there is
%   none. Call it before building a dialog and abort as if cancelled when
%   it returns something.
%
%   GUI.I_DLGREGISTER(parentfig, d) records D as the dialog now open on
%   PARENTFIG.
%
%   WHY THIS EXISTS. GUI.MYLISTDLG and GUI.MYTABLEDLG deliberately avoid
%   relying on WindowStyle='modal' for placement: on a multi-monitor setup
%   with mixed DPI, MATLAB re-centres modal dialogs onto the primary
%   monitor regardless of the Position given. What they use instead is
%   UIWAIT, and UIWAIT blocks only the CALLING code -- it keeps pumping
%   events, so the parent figure stays interactive. A second click on the
%   same toolbar button re-enters the callback and builds a second copy of
%   the dialog. That was visible in scgeatool on 'Show Cell States...' and
%   'Export/Save Data...', both of which reach GUI.MYLISTDLG.
%
%   The register is shared rather than one per dialog function because
%   GUI.MYLISTDLG delegates to GUI.MYTABLEDLG for very long lists. With
%   separate registers the delegating path would be unguarded: the caller
%   would find its own register empty, hand off, and open a second table
%   dialog.
%
%   No unregister call is needed. An entry holds the dialog handle, the
%   handle is deleted before the owning function returns, and dead entries
%   are pruned on the next call. If an error escapes with the dialog still
%   open the entry keeps holding, which is right -- the dialog is also
%   still on screen.
%
%   Only figures are tracked. A parentfig of [] registers nothing and is
%   never reported busy: with no parent there is no button to click twice.
%
% see also: GUI.MYLISTDLG, GUI.MYTABLEDLG

persistent openDlgs

if isempty(openDlgs)
    openDlgs = struct('parent', {}, 'dlg', {});
end
if ~isempty(openDlgs)
    openDlgs = openDlgs(arrayfun(@(e) pkg.i_isvalid(e.dlg), openDlgs));
end

existing = [];
if ~pkg.i_isvalid(parentfig)
    return;
end

if nargin < 2
    for k = 1:numel(openDlgs)
        if isequal(openDlgs(k).parent, parentfig)
            existing = openDlgs(k).dlg;
            i_raise(existing);
            return;
        end
    end
    return;
end

if pkg.i_isvalid(d)
    openDlgs(end+1) = struct('parent', parentfig, 'dlg', d);
end

end


function i_raise(d)
% Bring the open dialog forward, so a second click reads as "it is already
% up, here it is" rather than as a dead click.
try
    d.Visible = 'on';
    figure(d);
catch
    try
        focus(d);
    catch
        % neither is available; the dialog is on screen regardless
    end
end
end
