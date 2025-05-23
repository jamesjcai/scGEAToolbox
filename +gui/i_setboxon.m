function i_setboxon(~, ~, hFig)
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
    if ~isempty(axesHandles)
        
        if rand > 0.5
            paraset = {10, 0.5};
        else
            paraset = {15, 1.5};
        end
        
        for k = 1:length(axesHandles)
            ax = axesHandles(k);
            if strcmp(ax.Box, 'on')
                set(ax, 'Fontsize', paraset{1}, 'LineWidth', paraset{2});
                box(ax, 'off');
                %ax.XAxis.Visible="off";
                %ax.YAxis.Visible="off";
            else
                set(ax, 'Fontsize', paraset{1}, 'LineWidth', paraset{2});
                box(ax,'on');
                %ax.XAxis.Visible="on";
                %ax.YAxis.Visible="on";
            end
        end
    else
        gui.myHelpdlg(hFig, ['No figures available in the current ' ...
            'window. The ''box on/off'' command cannot be applied.']);
    end

end
