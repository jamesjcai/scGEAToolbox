function [answer] = myInputwin(~, ~, definput, parentfig)

answer = [];
hFig = uifigure("WindowStyle","modal",'Visible','off');

if ~isMATLABReleaseOlderThan('R2025a')
    try
        theme(hFig, parentfig.Theme.BaseColorStyle);
    catch ME
        disp(ME.message);
    end
end

hFig.Position(3)=0.75*hFig.Position(3);
hFig.Position(4)=0.75*hFig.Position(4);
dialogWidth = hFig.Position(3);
dialogHeight = hFig.Position(4);

% Compute center position relative to parentfig
parentPos = parentfig.Position;
dialogX = parentPos(1) + (parentPos(3) - dialogWidth) / 2;
dialogY = parentPos(2) + (parentPos(4) - dialogHeight) / 2;
hFig.Position(1)=dialogX;
hFig.Position(2)=dialogY;


g = uigridlayout(hFig,[3 3]);
g.RowHeight = {'fit','2x','fit'};
g.ColumnWidth = {'1x',75,75};

lbl = uilabel(g,"Text",'Paste List:');
lbl.Layout.Row = 1;
lbl.Layout.Column = 1;

txa = uitextarea(g);
txa.Layout.Row = 2;
txa.Layout.Column = [1 3];
txa.Value = compose(definput);

% txa.Position(4)=txa.Position(4)*2;
if nargout>0
    oktext = "OK";
else
    oktext = "Cancel";
end

btn = uibutton(g,"Text",oktext);
btn.Layout.Row = 3;
btn.Layout.Column = 2 + (nargout==0);
btn.ButtonPushedFcn = @(src,event) textEntered(src,event,btn);
%btn.Position(3) = 50;

if nargout > 0
    btn2 = uibutton(g,"Text","Cancel");
    btn2.Layout.Row = 3;
    btn2.Layout.Column = 3;
    btn2.ButtonPushedFcn = @(src,event) textEntered(src,event,btn2);
end

% gui.i_movegui2parent(hFig, parentfig);
drawnow;
pause(0.7);
hFig.Visible=true;
uiwait(hFig);

function textEntered(~,~,btn)
    if strcmp(btn.Text,oktext)
        y = true;
        answer = arrayfun(@(f) f.Value, txa, 'UniformOutput', false);
    else
        y = false;
        answer = {};
    end
    delete(hFig);
end

end

    % if nargin<1, parentfig = uifigure('Visible', 'off'); end
    % 
    % figure(parentfig);
    % % Example usage with multiple inputs
    % prompt = {'Enter matrix size:', 'Enter colormap name:'};
    % dlgtitle = 'Input';
    % % fieldsize = [1 45; 1 45]; % Not used in uifigure, but kept for consistency
    % definput = {'20', 'hsv'};


    %{
function answer = openInputDialog(prompt, dlgtitle, definput, parentfig)


    % Create modal dialog
    d = uifigure('Position', [dialogX, dialogY, dialogWidth, dialogHeight], ...
                 'Name', dlgtitle, 'WindowStyle', 'modal');

    % Store user input
    numFields = numel(prompt);
    fields = gobjects(numFields, 1);
    
    % Create UI elements
    for i = 1:numFields
        uilabel(d, 'Text', prompt{i}, ...
                   'Position', [20, dialogHeight - 40 - (i - 1) * 50, 310, 22]);
        fields(i) = uieditfield(d, 'text', ...
                                'Position', [20, dialogHeight - 65 - (i - 1) * 50, 310, 22], ...
                                'Value', definput{i});
    end

    % OK button
    okBtn = uibutton(d, 'Text', 'OK', ...
                     'Position', [100, 20, 60, 30], ...
                     'ButtonPushedFcn', @(btn, event) onOKButton(d, fields));

    % Cancel button
    cancelBtn = uibutton(d, 'Text', 'Cancel', ...
                         'Position', [190, 20, 60, 30], ...
                         'ButtonPushedFcn', @(btn, event) delete(d));

    % Wait for user input
    uiwait(d);
    
    % Retrieve data
    if isvalid(d) % If the dialog was not closed by user
        answer = arrayfun(@(f) f.Value, fields, 'UniformOutput', false);
        delete(d);
    else
        answer = {};
    end
end

function onOKButton(d, fields)
    % Resume UI execution when OK is pressed
    uiresume(d);
end
    %}