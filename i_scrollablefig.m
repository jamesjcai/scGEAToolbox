figure(1)
panel1 = uipanel('Parent',1);
panel2 = uipanel('Parent',panel1);
set(panel1,'Position',[0 0 .99 1]);
set(panel2,'Position',[-1 0 2 1]);
set(gca,'Parent',panel2);
h = image;    % <- Code for my own image goes here 
s = uicontrol('Style','Slider','Parent',1,...
    'Units','normalized','Position',[0 0 1 .025],...
    'Value',1,'Callback',{@slider_callback1,panel2});

function slider_callback1(src,~,arg1)
val = get(src,'Value');
set(arg1,'Position',[0 -val 1 2])
end