function [hFig] = sc_simpleabout(parentfig)

    % mfolder = fileparts(mfilename('fullpath'));
    % splashpng = '700813831-hero-1536x1536.png';
    % [im] = imread(fullfile(mfolder,'..','resources', splashpng));
    
mfolder = fileparts(mfilename('fullpath'));
if nargin<1
    parentfig=[];
end

    hFig = figure('MenuBar','none','ToolBar','none', ...
        'Name','About','NumberTitle','off','Color','k', ...
        'Position',[0 0 400 300],'Visible','off', ...
        'WindowStyle','modal', ...
        'DockControls','off', ...
        'Resize','off', ...
        'WindowButtonDownFcn',@(src, ~) close(src),...
        'WindowKeyPressFcn',@(src, ~) close(src));
    
    fa = axes('Parent',hFig,'Color','k', ...
        'XColor','k','YColor','k');
    %    'Units','Normalize','Position',[0 0 1 1]);    
    % set(gcf, 'WindowButtonDownFcn', @(src, event) close(gcf));
    %image(fa,im);
    [~, v1] = pkg.i_majvercheck;

    
    stfile =   fullfile(mfolder,'..','TIMESTAMP');
    if ~exist(stfile,'file')
        d = '';
    else
        fid=fopen(stfile,'r');
        d=textscan(fid,'%s');
        d=d{1}{1};
        fclose(fid);
    end

    % text(fa,0.75,0.2,'Loading...','Color','w');
    % plot(fa,[0.1 repmat(0.1, 1, r)],'-','LineWidth',5,'Color',[.7 .7 .7]);
    ylim(fa,[0 1]);
    xlim(fa,[1 10]);
    movegui(hFig,'center');
    set(fa,'Color','k','XColor','k','YColor','k');
    text(fa,0.2, 0.8,'SCGEATOOL','Color','w','FontSize',18);
    if ~isempty(d)
        text(fa,0.2,0.675,v1+" (build "+d+")",'Color',[.7 .7 .7],'FontSize',14);
    else
        text(fa,0.2,0.675,v1,'Color',[.7 .7 .7],'FontSize',14);
    end
    %text(fa,1.0,0.2,'James Cai (jcai@tamu.edu)','Color',[.7 .7 .7],'FontSize',14);
% try
%     if ~isempty(parentfig)
%         px_new = gui.i_getchildpos(parentfig,hFig);
%         movegui(hFig,px_new);
%     else
%         movegui(hFig,'center');
%     end
% catch
%     movegui(hFig, 'center');
% end
gui.i_movegui2parent(hFig, parentfig);

    hFig.Visible=true;
    
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

