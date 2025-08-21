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

    % Compute center position relative to parentfig
    if ~isempty(parentfig)
        parentPos = parentfig.Position;
        dialogX = parentPos(1) + (parentPos(3) - dialogWidth) / 2;
        dialogY = parentPos(2) + (parentPos(4) - dialogHeight) / 2;
        fig = uifigure('Position', ...
            [dialogX, dialogY, dialogWidth, dialogHeight], ...
            'WindowStyle','modal','visible','off', Name="", Icon="");
    else
        dialogX = 300;
        dialogY = 300;
        fig = uifigure('Position', ...
            [dialogX, dialogY, dialogWidth, dialogHeight], ...
            'WindowStyle','modal', 'visible','off', Name="", Icon="");
        movegui(fig, 'center');
    end

dimPanel = uipanel(fig, 'Position', [20 20 300 150], ...
    'Title',prompt);

if iscell(definput)
    definput = definput{1};
end
edit = uieditfield(dimPanel, 'Position',[20 70 260 22], 'Value', definput);

btnOk = uibutton(dimPanel, 'Text','OK', ...
                 'Position',[60 20 80 30], ...
                 'ButtonPushedFcn', @(~,~) uiresume(fig));
btnCancel = uibutton(dimPanel, 'Text','Cancel', ...
                     'Position',[160 20 80 30], ...
                     'ButtonPushedFcn', @(~,~) delete(fig));
% Show the dialog

if ~isMATLABReleaseOlderThan('R2025a')
    try
        theme(fig, parentfig.Theme.BaseColorStyle);
    catch
    end
end

fig.Visible = 'on';
uiwait(fig);

if isvalid(fig) % If the dialog was not closed by user
    answer = {edit.Value};
    delete(fig);
else
    answer = {};
end

% if strcmp(fig.SelectionType,'normal')
%     answer = edit.Value;
% else
%     answer = [];
% end


