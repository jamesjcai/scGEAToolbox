%hFigure=figure('Position',[100 100 500 300],'Visible','off');
hFigure=figure('Visible','off');

set(hFigure, 'MenuBar', 'none');
set(hFigure, 'ToolBar', 'none');

% set(f,'resizefcn',@myResizeFun);
movegui(hFigure,"center")
height = 25;
width = 120;
sz = hFigure.Position;
x = sz(3);
y = sz(4);
%x = mean( sz( [1, 3]));
%y = mean( sz( [2, 4]));
Position= [(x - width)/2, (y - height)/2, width, height];

tabgp = uitabgroup(hFigure,"Position",[.05 .05 .8 .7]);
tab1 = uitab(tabgp,"Title","Settings",'Tag',"1");
tab2 = uitab(tabgp,"Title","Options",'Tag',"2");
tabgp.SelectionChangedFcn=@displaySelection;



% % text(Position(1),Position(2),0,'Ready to explore.');
% a = uicontrol('style','push',...
%             'Parent',hFigure,...
%                  'position',Position,...
%                  'string','Import Data...',...
%                  'callback',{@gui.sc_openscedlg});
% % set(a,"Visible","off");
% Position(2)=Position(2)+24;
% %Position(3)=Position(3)*2;
% b = uicontrol('style','text',...
%             'Parent',hFigure,...
%             'FontSize',9,...
%             'position',Position,...
%             'string','Ready to explore.');


hAx = axes('Parent',hFigure,'Visible','on');
title(hAx, 'sce.title');
subtitle(hAx,'[genes x cells]');

hax1 = axes('Parent', tab1);
hax2 = axes('Parent', tab2);

h1=scatter(hax1,randn(300,1),rand(300,1));
h2=scatter(hax2,randn(300,1),rand(300,1),'rs');

%dt = datacursormode(hFigure);
%dt.UpdateFcn = {@i_myupdatefcnx};

drawnow;
set(hFigure,'Visible','on');
set(hAx,'Visible','on');



function displaySelection(src,event)
    t = event.NewValue;
    title = t.Title;
    disp("Viewing the " + title + " tab")
    % h1.Visible="on";
    % h2.Visible="off";
end


