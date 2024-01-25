function [name] = uiinputdlg(fig)
name = [];
height = 200;
width = 300;
sz = fig.Position;
x = sz(1) + sz(3)/2;
y = sz(2) + sz(4)/2;

%x = mean( sz( [1, 3]));
%y = mean( sz( [2, 4]));
%x = sz(1);
%y = sz(2);

f = uifigure(Visible="off",WindowStyle="modal");
f.Position= [x - width/2, y - height/2, width, height];
% hPrompt = uicontrol('Parent', f, 'Style', 'text', ...
%                    'String', 'Enter your name:');
% hNameInput = uicontrol('Parent', f, 'Style', 'edit', ...
%                       'Position', [50 50 150 25]);
% hOKButton = uicontrol('Parent', f, 'Style', 'pushbutton', ...
%                      'String', 'OK', ...
%                      'Callback', @handleInput);
b = uibutton(f,'ButtonPushedFcn', @handleInput, ...
    "Position",[10 10 80 20],'Text',"OK");

hNameInput = uilistbox(f);
hNameInput.Items = ["Red","Green","Blue"];

f.Visible='on';
%focus(f);
uiwait(f);

function handleInput(hObject, eventdata)
    %name = get(hNameInput, 'String');
    name = hNameInput.Value;
    delete(f);
end

end