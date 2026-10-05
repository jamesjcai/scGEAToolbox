function [fw] = myWaitbar(parentfig, fw, witherror, mesg, newmesg, f)
% MYWAITBAR - Open, update and close a progress dialog.
%
%   fw = gui.myWaitbar(parentfig)                      open
%   gui.myWaitbar(parentfig, fw, false, '', newmesg)   new message, spins
%   gui.myWaitbar(parentfig, fw, false, '', '', f)     fraction F done
%   gui.myWaitbar(parentfig, fw)                       close ("Finishing")
%   gui.myWaitbar(parentfig, fw, true)                 close after an error
%
% The dialog is always a UIPROGRESSDLG. It spins until a caller passes a
% fraction F, since without one nothing is known about how far along the work
% is; it used to be parked at 0.618, which read as stuck two thirds of the way.
% When PARENTFIG is not a uifigure the dialog goes on the open app window, or,
% with none open, on a small window of its own that closes along with it. The
% old WAITBAR figure could not spin at all.

if nargin < 6, f = []; end
if nargin < 5, newmesg = ''; end
if nargin < 4, mesg = ''; end
if nargin < 3 || isempty(witherror), witherror = false; end
if nargin < 1, parentfig = []; end

if nargin < 2 || isempty(fw)
    % The two defaults for MESG belong to different calls: 'Processing your
    % data...' when the dialog is created, 'Finishing' when it is closed.
    if isempty(mesg), mesg = 'Processing your data...'; end
    try
        [host, ownsHost] = i_hostfigure(parentfig);
        fw = uiprogressdlg(host, 'Title', 'Please wait...', ...
            'Message', mesg, 'Indeterminate', 'on');
    catch ME
        % No display (matlab -batch, a worker): run without a dialog. Every
        % later call is a no-op on the empty handle.
        warning('myWaitbar:noDialog', 'No progress dialog: %s', ME.message);
        fw = [];
        return;
    end
    if ownsHost
        addlistener(fw, 'ObjectBeingDestroyed', @(~, ~) delete(host));
    end
    fprintf('Processing your data...');
    fprintf('... ');
    tic;
    return;
end

if ~pkg.i_isvalid(fw) || ~isa(fw, 'matlab.ui.dialog.ProgressDialog')
    return;
end

if ~isempty(newmesg) && isempty(f)
    fw.Indeterminate = 'on';
    fw.Message = newmesg;
elseif isempty(newmesg) && ~isempty(f)
    fw.Indeterminate = 'off';
    fw.Value = f;
elseif ~isempty(newmesg) && ~isempty(f)
    fw.Indeterminate = 'off';
    fw.Message = newmesg;
    fw.Value = f;
else
    if ~witherror
        if isempty(mesg), mesg = 'Finishing'; end
        toc;
        fw.Indeterminate = 'off';
        fw.Value = 1;
        fw.Message = mesg;
        pause(1);
    end
    if pkg.i_isvalid(fw), close(fw); end
end
end

function [host, ownsHost] = i_hostfigure(parentfig)
% The uifigure to put the dialog on: PARENTFIG, the app window running the
% current callback, any open uifigure, or else a new small window, which
% OWNSHOST says to delete with the dialog. Another dialog's own window is
% never reused: it is deleted when that dialog closes.
ownsHost = false;
if gui.i_isuifig(parentfig)
    host = parentfig;
    return;
end
candidates = [gcbf; findall(groot, 'Type', 'figure', 'Visible', 'on')];
for k = 1:numel(candidates)
    if gui.i_isuifig(candidates(k)) && strcmp(candidates(k).Visible, 'on') && ...
            ~strcmp(candidates(k).Tag, 'myWaitbarHost')
        host = candidates(k);
        return;
    end
end
ownsHost = true;
host = uifigure('Name', 'Please wait...', 'Tag', 'myWaitbarHost', ...
    'Position', [0 0 420 140], 'Resize', 'off', 'Visible', 'off');
movegui(host, 'center');
if ~isempty(parentfig) && ishghandle(parentfig)
    newpos = gui.i_getchildpos(parentfig, host);
    if ~isempty(newpos), host.Position(1:2) = newpos; end
end
host.Visible = 'on';
end
