function [name] = uiinputdlg(fig,Items)

if nargin<2
    Items = ["aaa","bbb","ccc"];
end
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

b = uibutton(f,'ButtonPushedFcn', @handleInput, ...
    "Position",[10 10 80 20],'Text',"OK");

%hNameInput = uidropdown(f);
%hNameInput.Items = Items;


hNameInput = uilistbox(f);
hNameInput.Items = Items;
hNameInput.Multiselect = "on";

% ["Red","Green","Blue","Green","Blue","Green","Blue","Green","Blue","Green","Blue"];

f.Visible='on';
%focus(f);
uiwait(f);

function handleInput(hObject, eventdata)
    %name = get(hNameInput, 'String');
    name = hNameInput.Value;
    delete(f);
end

end