function i_stem3scatter(x,y,c,t)
if nargin<4, t=''; end
hFig=figure;
    stem3(x,y,c,'marker','none','color','m');
    hold on
    scatter3(x,y,zeros(size(y)),5,c,'filled');
    title(t);
    % hFig=gcf;
    tb = uitoolbar(hFig);
    gui.add_3dcamera(tb,'HgB_markers');
end
    
