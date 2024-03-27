function [fx] = sc_simpleabout2(parentfig)

    mfolder = fileparts(mfilename('fullpath'));
    splashpng = '700813831-hero-1536x1536.png';
    [im] = imread(fullfile(mfolder,'..','resources', splashpng));
    
if nargin<1
    parentfig=[];
end

    fx = figure('MenuBar','none','ToolBar','none', ...
        'Name','About','NumberTitle','off','Color','k', ...
        'Position',[0 0 400 300],'Visible','off', ...
        'WindowStyle','modal', ...
        'DockControls','off', ...
        'Resize','off', ...
        'WindowButtonDownFcn',@(src, ~) close(src),...
        'WindowKeyPressFcn',@(src, ~) close(src));
    
    fa = axes('Parent',fx,'Color','k', ...
        'XColor','k','YColor','k','Visible','off');
    %    'Units','Normalize','Position',[0 0 1 1]);    
    % set(gcf, 'WindowButtonDownFcn', @(src, event) close(gcf));
    %image(fa,im);

    %load earth.mat
    ih = image(im,'Parent',fa);
    %colormap(map)
    
imxpos = get(ih,'XData');
imypos = get(ih,'YData');
set(fa,'Unit','Normalized','Position',[0,0,1,1]);
figpos = get(fx,'Position');
figpos(3:4) = [imxpos(2) imypos(2)];
set(fx,'Position',figpos);

    % [~, v1] = pkg.i_majvercheck;
    v1 = pkg.i_getversionnum;

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
    % ylim(fa,[0 1]);
    % xlim(fa,[1 10]);
    movegui(fx,'center');
    set(fa,'Color','k','XColor','k','YColor','k');
    text(fa, 20, 50,'SCGEATOOL','Color','w','FontSize',16);
    if ~isempty(d)
        text(fa,20, 80, "Version "+v1+" (build "+d+")",'Color',[.7 .7 .7],'FontSize',12);
    else
        text(fa,20, 80,"Version "+v1,'Color',[.7 .7 .7],'FontSize',12);
    end
    %text(fa,1.0,0.2,'James Cai (jcai@tamu.edu)','Color',[.7 .7 .7],'FontSize',14);

    if ~isempty(parentfig)
        px_new = gui.i_getchildpos(parentfig,fx);
        movegui(fx,px_new);
    else
        movegui(fx,'center');
    end
    fx.Visible=true;
    
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

