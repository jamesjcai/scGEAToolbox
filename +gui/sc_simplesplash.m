function [fx] = sc_simplesplash(fx,r)

if nargin < 2, r = 1; end
if nargin < 1 || isempty(fx)
    fx = figure('MenuBar','none','ToolBar','none', ...
        'Name','','NumberTitle','off','Color','k', ...
        'Position',[0 0 400 300],'Visible','off','WindowStyle','normal');
    fa = axes('Parent',fx,'Color','k', 'XColor','k','YColor','k');
    % text(fa,0.75,0.2,'Loading...','Color','w');
    plot(fa,[0.1 repmat(0.1, 1, r)],'-','LineWidth',5,'Color',[.7 .7 .7]);
    ylim(fa,[0 1]);
    xlim(fa,[1 10]);
    movegui(fx,'center');
    set(fa,'Color','k','XColor','k','YColor','k');
    text(fa,0.2, 0.8,'SCGEATOOL','Color','w','FontSize',18);
    text(fa,0.2,0.675,'v24.3.3','Color',[.7 .7 .7],'FontSize',14);
    text(fa,1.0,0.2,'Loading...','Color',[.7 .7 .7],'FontSize',14);
    fx.Visible=true;
    hold on
else
    fa = findall(fx,'type','axes');    
    plot(fa,[0.1 repmat(0.1, 1, r)],'-','LineWidth',5,'Color',[.7 .7 .7]);    
end

% for k = 3:10
%     h = plot(fa, repmat(0.1, 1, k),'-','LineWidth',5,'Color',[.7 .7 .7]);    
%     drawnow;
%     pause(1)
% end
% 

% for k=10:-1:2
%     xlim(fa,[1 k]);
%     pause(.5);
%     drawnow;
% end

