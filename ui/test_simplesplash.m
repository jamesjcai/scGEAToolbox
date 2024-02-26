
fx = figure('MenuBar','none','ToolBar','none', ...
    'Name','Loading...','NumberTitle','off','Color','k', ...
    'Position',[0 0 350 250],'Visible','off','WindowStyle','modal');
fa = axes('Parent',fx,'Color','k', 'XColor','k','YColor','k');
% text(fa,0.75,0.2,'Loading...','Color','w');
h = plot(fa,[.1 .1],'-','LineWidth',5,'Color',[.7 .7 .7]);
hold on
ylim(fa,[0 1]);
xlim(fa,[1 10]);
movegui(fx,'center');
set(fa,'Color','k','XColor','k','YColor','k');
text(fa,0.2,0.8,'SCGEATOOL','Color','w');
fx.Visible=true;
for k = 3:10
    h = plot(fa, repmat(0.1, 1, k),'-','LineWidth',5,'Color',[.7 .7 .7]);    
    drawnow;
    pause(1)
end


% for k=10:-1:2
%     xlim(fa,[1 k]);
%     pause(.5);
%     drawnow;
% end

