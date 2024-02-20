f=figure('Position',[100 100 500 300]);
% set(f,'resizefcn',@myResizeFun);
movegui(f,"center")
height = 25;
width = 120;
sz = f.Position;
x = sz(3);
y = sz(4);
%x = mean( sz( [1, 3]));
%y = mean( sz( [2, 4]));
Position= [(x - width)/2, (y - height)/2, width, height];

% text(Position(1),Position(2),0,'Ready to explore.');
a = uicontrol('style','push',...
            'Parent',f,...
                 'position',Position,...
                 'string','Import Data...',...
                 'callback',{@gui.sc_openscedlg});
% set(a,"Visible","off");
Position(2)=Position(2)+24;
%Position(3)=Position(3)*2;
b = uicontrol('style','text',...
            'Parent',f,...
            'FontSize',9,...
            'position',Position,...
            'string','Ready to explore.');



hAx = axes('Parent', f,'Visible','off');
title(hAx, 'sce.title');
subtitle(hAx,'[genes x cells]');
dt = datacursormode(f);
dt.UpdateFcn = {@i_myupdatefcnx};
return;


%%

% fig = figure;
fig=figure('Position',[300 300 500 300]);

fig_pos = get(fig, 'Position'); % [left bottom width height]
fig_width = fig_pos(3);
fig_height = fig_pos(4);

btn_width = 100; % Adjust as needed
btn_height = 25; % Adjust as needed
btn_x = (fig_width - btn_width) / 2;
btn_y = (fig_height - btn_height) / 2;

button = uicontrol(...
    'Style', 'pushbutton',...
    'Units', 'pixels',...
    'Position', [btn_x btn_y btn_width btn_height],...
    'String', 'Click Me!',... % Customize label
    'FontSize', 14,... % Customize font size
    'BackgroundColor', 'blue',... % Customize background color
    'ForegroundColor', 'white',... % Customize text color
    'Callback', {@gui.sc_openscedlg}); % Assign callback function


% uibutton(fig)
set(fig,'resizefcn',{@myResizeFun,button});
plot(rand(300,1));


%%
fig = figure('Position',[561 497 333 239]);
h = uihtml(fig);
h.HTMLSource = '<p style="font-family:arial;"><b><span style="color:red;">Hello</span> <u>World</u>!</b></p>';

%%
% fig = figure;
% im = uiimage(fig);
% im.ImageSource = "peppers.png";

 f= figure()
 a = axes('Position',[0 0 0.1 0.1],'Units','Normalized');
 imshow('peppers.png','Parent',a);
 uicontrol('Style','text','Units', ...
     'Normalized','Position', ...
     [0.1 0.1 0.3 0.6],'String','Example');

%%


function myResizeFun(src,events,butt)
src
events
disp('aaa')

fig_pos = get(src, 'Position'); % [left bottom width height]
fig_width = fig_pos(3);
fig_height = fig_pos(4);

btn_width = 100; % Adjust as needed
btn_height = 25; % Adjust as needed
btn_x = (fig_width - btn_width) / 2;
btn_y = (fig_height - btn_height) / 2;

set(butt,'Position',[btn_x btn_y btn_width btn_height]);
end