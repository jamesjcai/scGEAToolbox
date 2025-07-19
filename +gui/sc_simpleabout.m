function [hFig] = sc_simpleabout(parentfig, im2)
    if nargin < 2 || isempty(im2)
        [~, ~, ~, im2] = pkg.i_majvercheck;
    end

    mfolder = fileparts(mfilename('fullpath'));


    splashdir = fullfile(mfolder, '..','assets', 'Images', 'splash_folder');
    a = dir(splashdir);
    idx = 2+randi(length(a)-2);
    pngfilename = a(idx).name;
    splashpng = fullfile(mfolder, '..','assets', 'Images','splash_folder', pngfilename);

    % splashpng = 'splash.png';
    % fullfile(mfolder,'..','assets', 'Images', splashpng)
    % splashpng = 'thumbnail_IMG_3621.jpg';
    [im] = imread(splashpng);
    if ~isempty(im2) && license('test','image_toolbox') && ~isempty(which('imfuse')) 
        im = imfuse(im, im2, 'blend');
    end
    if nargin<1, parentfig=[]; end

    
    iconfile = fullfile(mfolder, '..','assets', 'Images', 'icon_16.png');

    done = false;
    try
        hFig = figure('MenuBar','none','ToolBar','none', ...
            'Name','About','NumberTitle','off','Color','k', ...
            'Position',[0 0 400 300],'Visible','off', ...
            'WindowStyle','modal', ...
            'DockControls','off', ...
            'Resize','off', ...
            'WindowButtonDownFcn',@(src, ~) close(src),...
            'WindowKeyPressFcn',@(src, ~) close(src), 'Icon', iconfile);
        done = true;
    catch ME
        % disp(ME.message)
    end

    if ~done
        hFig = figure('MenuBar','none','ToolBar','none', ...
            'Name','About','NumberTitle','off','Color','k', ...
            'Position',[0 0 400 300],'Visible','off', ...
            'WindowStyle','modal', ...
            'DockControls','off', ...
            'Resize','off', ...
            'WindowButtonDownFcn',@(src, ~) close(src),...
            'WindowKeyPressFcn',@(src, ~) close(src));
    end

    fa = axes('Parent', hFig, 'Color','k', ...
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
    figpos = get(hFig,'Position');
    figpos(3:4) = [imxpos(2) imypos(2)];
    set(hFig,'Position',figpos);

    % [~, v1] = pkg.i_majvercheck;
    v1 = pkg.i_get_versionnum;

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
    movegui(hFig,'center');
    set(fa,'Color','k','XColor','k','YColor','k');
    % text(fa, 40, 10,'since 2024/10', 'Color', [.5 .5 .5]);
    text(fa, 20, 50,'SCGEATOOL','Color','w','FontSize',16);
    if ~isempty(d)
        text(fa, 20, 80, "Version "+v1+" (build "+d+")",'Color',[.7 .7 .7],'FontSize',12);
    else
        text(fa, 20, 80,"Version "+v1,'Color',[.7 .7 .7],'FontSize', 12);
    end
    %text(fa,1.0,0.2,'James Cai (jcai@tamu.edu)','Color',[.7 .7 .7],'FontSize',14);


% text(fa, 20, 180, sprintf("DISCLAIMER:\n\n" + ...
%     "The software is distributed ""as is,"" without warranty of any kind,\n" + ...
%     "either express or implied, including but not limited to the \n" + ...
%     "warranties of merchantability, fitness for a particular purpose, and\n" + ...
%     "non-infringement. In no event shall the authors or Texas A&M \n" + ...
%     "University be liable for any claim, damages, or other liability\n" + ...
%     "arising from the use of this software."), ...
%     'Color',[.7 .7 .7], 'FontSize', 9);

disclaimer_text = 'DISCLAIMER: The software is distributed "as is," without warranty of any kind, either express or implied, including but not limited to the warranties of merchantability, fitness for a particular purpose, and non-infringement. In no event shall the authors or Texas A&M University be liable for any claim, damages, or other liability arising from the use of this software.';

annotation('textbox', [0.05 0.1 0.9 0.5], ...
    'String', disclaimer_text, ...
    'FontSize', 10, ...    
    'Color',[.7 .7 .7],...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle','EdgeColor','None');

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

