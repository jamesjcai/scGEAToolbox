function sc_uitabgrpfig(sce, glist, parentfig)

if nargin<3, parentfig = []; end
if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
    p = parentfig.Position;
    cx = [p(1)+p(3)/2 p(2)+p(4)/2];
end

import mlreportgen.ppt.*;

pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, 'resources', 'myTemplate.pptx');


hFig=figure("Visible","off");
hFig.Position(3) = hFig.Position(3) * 1.8;

if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
    px = hFig.Position;
    px_new = [cx(1)-px(3)/2 cx(2)-px(4)/2];
    
    % if px_new(1)<0
    %     ss = get(0, 'screensize');
    %     px_new(1)=px_new(1)-ss(4); %ss(3); 
    % end
else
    px_new = [];
end

n = length(glist);
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

tabgp = uitabgroup();

idx = 1;
focalg = glist(idx);

for k=1:n
    c = sce.X(sce.g == glist(k), :);
    if issparse(c), c = full(c); end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',glist(k)));
    
    %{
    t = tiledlayout(1,2,'Parent',tab{k});
    ax1 = nexttile;
    hpl{k,1} = scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled','Parent', ax1);
    ax2 = nexttile;
    hpl{k,2} = scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled','Parent', ax2);
    %}
    
    ax0{k} = axes('parent',tab{k});
    ax{k,1} = subplot(1,2,1);
    scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled');
    ax{k,2} = subplot(1,2,2);
    scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled');
    stem3(sce.s(:,1), sce.s(:,2), c, 'marker', 'none', 'color', 'm');
    hold on;
    scatter3(sce.s(:,1), sce.s(:,2), zeros(size(sce.s(:,2))), 5, c, 'filled');

    title(ax{k,1}, glist(k));
    subtitle(ax{k,1}, gui.i_getsubtitle(c));
    title(ax{k,2}, glist(k));
    subtitle(ax{k,2}, gui.i_getsubtitle(c));

    %{
    ax = axes('Parent', tab{k});
    hold(ax, 'on');
    ax1 = subplot(1,2,1,tab{k});
    hold(ax1,'on')
    hpl{k,1} = scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled','Parent', ax1);
    ax2 = subplot(1,2,2,tab{k});
    hold(ax2,'on')
    hpl{k,2} = scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled','Parent', ax2);
    %}

    %{
    hax{k} = axes('Parent', tab{k});
    if size(sce.s,2)>=3
        hpl{k} = scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled','Parent', hax{k});
    else
        hpl{k} = scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled','Parent', hax{k});
    end
    title(hax{k}, targetg(k));
    subtitle(hax{k}, gui.i_getsubtitle(c));
    %}


    gui.i_setautumncolor(c, a, true, any(c==0));
end
  
tabgp.SelectionChangedFcn=@displaySelection;

tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
% b=allchild(tb0)
% tb = uitoolbar(hFig);
% copyobj(b(4),tb);
% delete(tb0);

uipushtool(tb, 'Separator', 'off');
% pkg.i_addbutton2fig(tb, 'off', [], "IMG00107.GIF", " ");
% pkg.i_addbutton2fig(tb, 'off', @i_linksubplots, 'plottypectl-rlocusplot.gif', 'Link subplots');
pkg.i_addbutton2fig(tb, 'on',  @i_genecards, 'fvtool_fdalinkbutton.gif', 'GeneCards...');
pkg.i_addbutton2fig(tb, 'on', {@i_PickColorMap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
%pkg.i_addbutton2fig(tb, 'off', @i_RescaleExpr, 'IMG00074.GIF', 'Rescale expression level [log2(x+1)]');
%pkg.i_addbutton2fig(tb, 'off', @i_ResetExpr, 'plotpicker-geobubble2.gif', 'Reset expression level');
% pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'off', @i_savemainfig, "powerpoint.gif", 'Save Figure to PowerPoint File...');

%gui.add_3dcamera(tb);
if isempty(px_new)
    movegui(hFig,'center');
else    
    movegui(hFig, px_new);
end
drawnow;
hFig.Visible=true;

    function i_savemainfig(~,~)
        answer = questdlg('Export to PowerPoint?');
        if ~strcmp(answer,'Yes'), return; end

        fw=gui.gui_waitbar_adv;
            OUTppt = [tempname, '.pptx'];
            ppt = Presentation(OUTppt, pth);
            open(ppt);
            images=cell(n,1);
            warning off
        for kx=1:n
            gui.gui_waitbar_adv(fw,kx./n,"Processing "+glist(kx)+" ...");
            images{kx} = [tempname, '.png'];
            tabgp.SelectedTab=tab{kx};
            saveas(tab{kx},images{kx});
            slide3 = add(ppt, 'Small Title and Content');
            replace(slide3, 'Title', glist(kx));
            replace(slide3, 'Content', Picture(images{kx}));        
        end
            close(ppt);
            rptview(ppt);      
            gui.gui_waitbar_adv(fw);
    end

    % function i_linksubplots(~,~)        
    %     hlink = linkprop([ax{idx,1},ax{idx,2}],{'CameraPosition','CameraUpVector'});
    % end

    function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        % disp("Viewing gene " + txt);
        [~,idx]=ismember(txt,glist);
        focalg = glist(idx);
    end

    function i_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end

end



    function i_PickColorMap(~, ~, c)
        list = {'parula', 'turbo', 'hsv', 'hot', 'cool', 'spring', ...
            'summer', 'autumn (default)', ...
            'winter', 'jet'};
        [indx, tf] = listdlg('ListString', list, 'SelectionMode', 'single', ...
            'PromptString', 'Select a colormap:');
        if tf == 1
            a = list{indx};
            if strcmp(a, 'autumn (default)')
                a = 'autumn';
            end
            gui.i_setautumncolor(c, a);
            setpref('scgeatoolbox', 'prefcolormapname', a);
        end
    end
