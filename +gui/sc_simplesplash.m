function [fx, v1] = sc_simplesplash(fx, r)
    if nargin<1, fx = []; end
    v1 = '';
    % https://www.mathworks.com/matlabcentral/answers/92259-how-do-i-make-a-splash-screen-for-my-matlab-gui-application
    mfolder = fileparts(mfilename('fullpath'));
    % splashpng = 'splash.png';
    splashdir = fullfile(mfolder, '..','assets', 'Images', 'splash_folder');
    a = dir(splashdir);

    if length(a)<=2, return; end

    idx = 2+randi(length(a)-2);
    pngfilename = a(idx).name;
    splashpng = fullfile(mfolder, '..','assets', 'Images','splash_folder', pngfilename);
   
    if ~isfile(splashpng)
        error('Splash image file not found: %s', splashpng);
    end    
    [im] = imread(splashpng);
    

if nargin < 2, r = 0.1; end
% if r>1, r=r/10; end
if nargin < 1 || isempty(fx)
    fx = figure('MenuBar','none','ToolBar','none', ...
        'Name','','NumberTitle','off','Color','k', ...
        'Position',[0 0 400 300],'Visible','off', ...
        'WindowStyle','normal', ...
        'DockControls','off', ...
        'Resize','off');
    
    % fa = axes('Parent',fx,'Color','k', ...
    %     'XColor','k','YColor','k','Visible','off');
  
        fa = axes('Parent',fx, 'Visible','off');
        ih = image(im,'Parent',fa);
    %colormap(map)
    
    imxpos = get(ih,'XData');
    imypos = get(ih,'YData');
    set(fa,'Unit','Normalized','Position',[0,0,1,1]);
    figpos = get(fx,'Position');
    figpos(3:4) = [imxpos(2) imypos(2)];
    set(fx,'Position',figpos);


    % text(fa,0.75,0.2,'Loading...','Color','w');
    %plot(fa,[255 repmat(255, 1, r)],'-','LineWidth',5,'Color',[.7 .7 .7]);
    hold on
    r=0.1;
    X=10:390;
    x=X(1:round(length(X)*r));
    y=290*ones(size(x));
    plot(fa,x,y,'-','LineWidth',4,'Color',[.7 .7 .7]);

    %ylim(fa,[0 1]);
    %xlim(fa,[1 10]);
    movegui(fx,'center');
    % set(fa,'Color','k','XColor','k','YColor','k');
    text(fa,20, 50,'SCGEATOOL','Color','w','FontSize',16);

    v1  = pkg.i_get_versionnum;
    text(fa,20, 80, v1, 'Color',[.7 .7 .7],'FontSize',12);
    text(fa,20, 270,'Loading...','Color',[.7 .7 .7],'FontSize',12);
    fx.Visible=true;
    hold on
else
    fa = findall(fx,'type','axes');    
    X=10:390;
    x=X(1:round(length(X)*r));
    y=290*ones(size(x));
    plot(fa,x,y,'-','LineWidth',4,'Color',[.7 .7 .7]);
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

