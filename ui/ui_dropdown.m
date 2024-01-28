function [name] = ui_dropdown(Parentfig,Items)

if nargin<2
    Items = ["aaa","bbb","ccc"];
end
name = [];
height = 200;
width = 300;
sz = Parentfig.Position;
x = sz(1) + sz(3)/2;
y = sz(2) + sz(4)/2;

Fig = uifigure(Visible="off",WindowStyle="modal");
Fig.Position= [x - width/2, y - height/2, width, height];

b = uibutton(Fig,'ButtonPushedFcn', @handleInput, ...
    "Position",[10 10 80 20],'Text',"OK");

hNameInput = uidropdown(Fig);
hNameInput.Items = Items;

%hNameInput = uilistbox(Fig);
%hNameInput.Items = Items;
%hNameInput.Multiselect = "on";

% ["Red","Green","Blue","Green","Blue","Green","Blue","Green","Blue","Green","Blue"];

Fig.Visible='on';
uiwait(Fig);

function handleInput(hObject, eventdata)
    %name = get(hNameInput, 'String');
    name = hNameInput.Value;
    delete(Fig);
end

end