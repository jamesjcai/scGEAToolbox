function [answer] = i_inputdlg(prompt, definput, parentfig)

if nargin<1, prompt = 'Enter Info'; end
if nargin<2, definput = ''; end
if nargin<3, parentfig = []; end

%{
f=uifigure;
gui.i_inputdlg('prompt', 'dlgtitle', 'definput', f)
%}

     dialogWidth = 340;
     dialogHeight = 100 + 90; % Adjust height based on the number of inputs

    pos = gui.i_centerdlgpos(parentfig, [dialogWidth, dialogHeight]);
    % WindowStyle='modal' is intentionally omitted: deleting a modal uifigure
    % hands focus back to the MATLAB desktop rather than to parentfig, which
    % drops the caller's window behind other windows and makes it look as if
    % the main figure vanished.  uiwait(fig) below still blocks the calling
    % code, so the dialog is functionally modal.  Same reasoning as in
    % gui.myListdlg.
    fig = uifigure('Position', pos, ...
        'Visible', 'off', Name="", Icon="");
    % Route the window's X button through cancelFcn so the dialog is never
    % closed by the default closereq, which falls back to close('force') --
    % and closes every open figure -- whenever gcbf is empty.
    fig.CloseRequestFcn = @(~,~) cancelFcn();

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

drawnow;
pause(0.7);
fig.Visible = 'on';
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
