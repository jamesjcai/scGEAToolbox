function i_changefontsize(~, ~, hFig)
    %{
        ax = get(gca, 'FontSize') - 1;
        if ax <= 5, ax = 15; end
        set(gca, 'FontSize', ax);
        %      ax=get(gca,'LabelFontSizeMultiplier')/1.1;
        %      if ax<=0.5, ax=3; end
        %      set(gca,'LabelFontSizeMultiplier')
    %}
    if nargin < 3, hFig = gcf; end
    axesHandles = findall(hFig, 'Type', 'axes');
    for k = 1:length(axesHandles)
        ax = axesHandles(k);
        fz = get(ax, 'FontSize') - 1;
        if fz <= 5, fz = 15; end
        set(ax, 'FontSize', fz);
    end

end
