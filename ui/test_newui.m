f=figure('Position',[10 10 200 300]);
movegui(f,"center")
height = 25;
width = 120;
sz = f.Position;
x = sz(3);
y = sz(4);
%x = mean( sz( [1, 3]));
%y = mean( sz( [2, 4]));
Position= [(x - width)/2, (y - height)/2, width, height];


a = uicontrol('style','push',...
            'Parent',f,...
                 'position',Position,...
                 'string','Import Data...',...
                 'callback',{@gui.sc_openscedlg});



%%

% fig = figure;
fig=figure('Position',[10 10 200 300]);

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
