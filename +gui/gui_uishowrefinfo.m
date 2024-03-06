function [txt,T] = gui_uishowrefinfo(reftarget)

txt = [];
pw1 = fileparts(mfilename('fullpath'));
fname = fullfile(pw1, '..','resources','refinfo.txt');
fid=fopen(fname,'r');
T=textscan(fid,'%s%s','Delimiter','\t');
fclose(fid);
reftag=string(T{:,1});

idx=find(reftarget==reftag);
if ~isempty(idx)
    txt=T{:,2}{idx};
end

% if ~isempty(txt)
%     fprintf('%s\n%s\n', reftarget, txt);
%     uiwait(helpdlg(txt, reftarget));
% end

fig = uifigure;
g = uigridlayout(fig,[3 1]);
g.RowHeight = {'fit','1x','fit'};
g.ColumnWidth = {'1x'};

lbl = uilabel(g,"Text",reftarget);
lbl.Layout.Row = 1;
lbl.Layout.Column = 1;
txa = uitextarea(g);
txa.Layout.Row = 2;
txa.Layout.Column = 1;
txa.Value = txt;
% txa.Position(4)=txa.Position(4)*2;
btn = uibutton(g,"Text","Submit","Enable","off");
btn.Layout.Row = 3;
btn.Layout.Column = 1;

txa.ValueChangedFcn = @(src,event) textEntered(src,event,btn);
end

function textEntered(src,event,btn)
val = src.Value;
btn.Enable = "off";
% Check each element of text area cell array for text
for k = 1:length(val)
    if ~isempty(val{k})
        btn.Enable = "on";
        break
    end
end
end

