function sc_uitabgrpfig(sce, targetg)

hFig=figure("Visible","off");       

n=length(targetg);
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

tabgp = uitabgroup();
focalg = targetg(1);
for k=1:n
    %tab{k} = uitab(group, 'Title', sprintf('Tab%d',k));
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',targetg(k)));
    hax{k} = axes('Parent', tab{k});
    c = sce.X(sce.g == targetg(k), :);
    if issparse(c), c = full(c); end
    if size(sce.s,2)>=3
        hpl{k} = scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled','Parent', hax{k});
    else
        hpl{k} = scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled','Parent', hax{k});
    end
    title(hax{k}, targetg(k));
    subtitle(hax{k}, gui.i_getsubtitle(c));
    gui.i_setautumncolor(c, a, true, any(c==0));
end
tabgp.SelectionChangedFcn=@displaySelection;

tb = uitoolbar(hFig);
%pkg.i_addbutton2fig(tb, 'off', @gui.i_linksubplots, 'plottypectl-rlocusplot.gif', 'Link subplots');
pkg.i_addbutton2fig(tb, 'on',  @i_genecards, 'fvtool_fdalinkbutton.gif', 'GeneCards...');
pkg.i_addbutton2fig(tb, 'on', {@i_PickColorMap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
%pkg.i_addbutton2fig(tb, 'off', @i_RescaleExpr, 'IMG00074.GIF', 'Rescale expression level [log2(x+1)]');
%pkg.i_addbutton2fig(tb, 'off', @i_ResetExpr, 'plotpicker-geobubble2.gif', 'Reset expression level');
pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
%gui.add_3dcamera(tb);
drawnow;
hFig.Visible=true;


    function displaySelection(src,event)
        t = event.NewValue;
        txt = t.Title;
        disp("Viewing the " + txt + " tab");
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
