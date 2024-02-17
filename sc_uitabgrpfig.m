function sc_uitabgrpfig(sce, targetg)

hFig=figure("Visible","off");       
hFig.Position(3) = hFig.Position(3) * 1.8;
n=length(targetg);
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

tabgp = uitabgroup();

idx = 1;
focalg = targetg(idx);

for k=1:n
    c = sce.X(sce.g == targetg(k), :);
    if issparse(c), c = full(c); end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',targetg(k)));
    
    %{
    t = tiledlayout(1,2,'Parent',tab{k});
    ax1 = nexttile;
    hpl{k,1} = scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled','Parent', ax1);
    ax2 = nexttile;
    hpl{k,2} = scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled','Parent', ax2);
    %}
    
    axes('parent',tab{k});
    ax{k,1} = subplot(1,2,1);
    scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled');
    ax{k,2} = subplot(1,2,2);
    scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled');
    stem3(sce.s(:,1), sce.s(:,2), c, 'marker', 'none', 'color', 'm');
    hold on;
    scatter3(sce.s(:,1), sce.s(:,2), zeros(size(sce.s(:,2))), 5, c, 'filled');

    title(ax{k,1}, targetg(k));
    subtitle(ax{k,1}, gui.i_getsubtitle(c));
    title(ax{k,2}, targetg(k));
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
pkg.i_addbutton2fig(tb, 'off', [], "IMG00107.GIF", " ");
% tb = uitoolbar(hFig);
% pkg.i_addbutton2fig(tb, 'off', @i_linksubplots, 'plottypectl-rlocusplot.gif', 'Link subplots');
pkg.i_addbutton2fig(tb, 'on',  @i_genecards, 'fvtool_fdalinkbutton.gif', 'GeneCards...');
pkg.i_addbutton2fig(tb, 'on', {@i_PickColorMap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
%pkg.i_addbutton2fig(tb, 'off', @i_RescaleExpr, 'IMG00074.GIF', 'Rescale expression level [log2(x+1)]');
%pkg.i_addbutton2fig(tb, 'off', @i_ResetExpr, 'plotpicker-geobubble2.gif', 'Reset expression level');
pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
%gui.add_3dcamera(tb);
movegui(hFig,'center');
drawnow;
hFig.Visible=true;


    function i_linksubplots(~,~)        
        hlink = linkprop([ax{idx,1},ax{idx,2}],{'CameraPosition','CameraUpVector'});
    end

    function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        % disp("Viewing gene " + txt);
        [~,idx]=ismember(txt,targetg);
        focalg = targetg(idx);
    end

    function i_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg));
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
