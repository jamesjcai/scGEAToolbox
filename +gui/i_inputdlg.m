function [answer] = i_inputdlg(prompt, definput, parentfig, modal)
%I_INPUTDLG Ask for one value, in a dialog centered on parentfig.
%
%   answer = gui.i_inputdlg(prompt, definput, parentfig, modal)
%
%   MODAL (default true) blocks input to the parent while the dialog is up.
%   UIWAIT alone blocks only the calling code, so without it the parent
%   stays clickable and a second click re-enters the callback. Pass false
%   for a dialog that must leave the parent usable.
%
%   ANSWER is a 1-by-1 cell holding the text, or {} when cancelled.
%
%   See also GUI.MYLISTDLG, GUI.I_RAISEFIG, GUI.I_CENTERDLGPOS.

if nargin<1, prompt = 'Enter Info'; end
if nargin<2, definput = ''; end
if nargin<3, parentfig = []; end
if nargin<4 || isempty(modal), modal = true; end

% At most one of these per parent. WINDOWSTYLE='modal' is not set until the
% dialog is on screen (see below), and a caller passing modal=false never
% gets it at all, so there is always a gap in which the parent -- and the
% toolbar button that opened this -- stays clickable. GUI.I_DLGREGISTER
% raises the dialog already up; returning {} then makes callers abort on
% their usual cancel branch.
if ~isempty(gui.i_dlgregister(parentfig))
    answer = {};
    return;
end

%{
f=uifigure;
gui.i_inputdlg('prompt', 'dlgtitle', 'definput', f)
%}

     dialogWidth = 340;
     dialogHeight = 100 + 90; % Adjust height based on the number of inputs

    pos = gui.i_centerdlgpos(parentfig, [dialogWidth, dialogHeight]);
    % Created non-modal and switched over once it is on screen, further
    % down. Two things used to keep WindowStyle='modal' out of here, and
    % both are now handled: deleting a modal uifigure hands focus to the
    % MATLAB desktop rather than to parentfig, which GUI.I_RAISEFIG at the
    % end of this function undoes; and making a figure modal before it is
    % realized lets MATLAB re-centre it onto the primary monitor, which is
    % why the switch happens after it is shown, as GUI.MYTABLEDLG does.
    fig = uifigure('Position', pos, ...
        'Visible', 'off', Name="", Icon="");
    % Route the window's X button through cancelFcn so the dialog is never
    % closed by the default closereq, which falls back to close('force') --
    % and closes every open figure -- whenever gcbf is empty.
    fig.CloseRequestFcn = @(~,~) cancelFcn();

gui.i_dlgregister(parentfig, fig);

if ~isMATLABReleaseOlderThan("R2025a")
        try
            % fig.Theme.BaseColorStyle = parentfig.Theme.BaseColorStyle;
            theme(fig, parentfig.Theme.BaseColorStyle);
        catch ME
            disp(ME.message);
        end
    end

dimPanel = uipanel(fig, 'Position', [20 20 300 150], ...
'Title', prompt);

if iscell(definput)
    definput = definput{1};
end
edit = uieditfield(dimPanel, 'Position', [20 70 260 22], 'Value', definput);

btnOk = uibutton(dimPanel, 'Text','OK', ...
                 'Position',[60 20 80 30], ...
                 'ButtonPushedFcn', @(~,~) uiresume(fig));
btnCancel = uibutton(dimPanel, 'Text','Cancel', ...
                     'Position',[160 20 80 30], ...
                     'ButtonPushedFcn', @(~,~) cancelFcn());
% Show the dialog

if ~isMATLABReleaseOlderThan('R2025a')
    try
        theme(fig, parentfig.Theme.BaseColorStyle);
    catch
        % theme() may not exist or parent has no Theme property; skip styling
    end
end

% No blind PAUSE before this point. PAUSE pumps the event queue while the
% dialog is still invisible and still non-modal, which is precisely when a
% second click on the opening button used to slip through and build a
% duplicate. DRAWNOW alone realizes the components, and the register above
% covers what is left. See GUI.MYTABLEDLG for the same sequence.
drawnow;
fig.Visible = 'on';

% Setting Visible, setting WindowStyle and DRAWNOW each pump the event
% queue with the dialog already on screen, so Cancel or the window's X can
% delete FIG between any two of the steps below. Check before every touch
% rather than assume it survived the previous one.
if modal && pkg.i_isvalid(fig)
    fig.WindowStyle = 'modal';
    if pkg.i_isvalid(fig), drawnow; end
    % Belt and braces against the re-centering described above; a no-op
    % when nothing moved it.
    if pkg.i_isvalid(fig) && ~isequal(round(fig.Position), round(pos))
        fig.Position = pos;
    end
end
if ~pkg.i_isvalid(fig)
    % Dismissed before it finished coming up; same result as Cancel.
    answer = {};
    gui.i_raisefig(parentfig);
    return;
end
focus(edit);
uiwait(fig);

if pkg.i_isvalid(fig) % If the dialog was not closed by user
    answer = {edit.Value};
    uiresume(fig);
    delete(fig);
else
    answer = {};
end

% Bring the caller's window back up: closing a dialog leaves the focus with
% whichever window Windows picks next, which is often not parentfig.
gui.i_raisefig(parentfig);

function cancelFcn()
        uiresume(fig);
        delete(fig);
    end
end

% if strcmp(fig.SelectionType,'normal')
%     answer = edit.Value;
% else
%     answer = [];
% end
