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

gui.i_movegui2parent(hFig, parentfig);


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
% btn.Position(3) = 50;

if nargout > 0
    btn2 = uibutton(g,"Text","Cancel");
    btn2.Layout.Row = 3;
    btn2.Layout.Column = 3;
    btn2.ButtonPushedFcn = @(src,event) textEntered(src,event,btn2);
end

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
