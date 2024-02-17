function [f]=sc_uitabgrpfig(sce, targetg)

f=figure;
n=length(targetg);
a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

group = uitabgroup();

    for k=1:n
        %tab{k} = uitab(group, 'Title', sprintf('Tab%d',k));
        tab{k} = uitab(group, 'Title', sprintf('%s',targetg(k)));
        hax{k} = axes('Parent', tab{k});
        c = sce.X(sce.g == targetg(k), :);
        if issparse(c), c = full(c); end
        if size(sce.s,2)>=3
            hpl{k} = scatter3(sce.s(:,1), sce.s(:,2), sce.s(:,3), 5, c, 'filled','Parent', hax{k});
        else
            hpl{k} = scatter(sce.s(:,1), sce.s(:,2), 5, c, 'filled','Parent', hax{k});
        end
        title(hax{k}, targetg(k))
        subtitle(hax{k}, gui.i_getsubtitle(c));
        gui.i_setautumncolor(c, a, true, any(c==0));    
    end

end

