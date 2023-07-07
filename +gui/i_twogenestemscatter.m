function i_twogenestemscatter(sce,g1,g2)

if ismember(upper(g1),upper(sce.g)) && ismember(upper(g2),upper(sce.g)) 

        hFig=figure('Visible','off');
        subplot(1,2,1);
            sc_scattermarker(sce.X,upper(sce.g),sce.s, ...
                upper(g1),1,[],false);
            box on
        subplot(1,2,2);
            sc_scattermarker(sce.X,upper(sce.g),sce.s, ...
                upper(g2),1,[],false);
            box on
        hFig.Position(3) = hFig.Position(3) * 1.55;
        set(hFig,'visible','on');
%         evalin('base', 'h=findobj(gcf,''type'',''axes'');');
%         evalin('base', 'hlink = linkprop(h(1:2),{''CameraPosition'',''CameraUpVector''});');
%         evalin('base', 'rotate3d on');
else
    disp('aaa');
end