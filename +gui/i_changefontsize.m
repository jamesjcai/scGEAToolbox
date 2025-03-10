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
    if ~isempty(axesHandles)
        hasax = false;

        for k = 1:length(axesHandles)
            ax = axesHandles(k);
            if ~isempty(ax)
                fz = get(ax, 'FontSize') - 1;
                if fz <= 5, fz = 15; end
                set(ax, 'FontSize', fz);
                hasax = true;
            end
        end
        if ~hasax
            gui.myHelpdlg(hFig, 'No plots available in the current window. Unable to change the font size', '');
        end
    else
        gui.myHelpdlg(hFig, 'No figures available in the current window. Unable to change the font size.', '');
    end
end
