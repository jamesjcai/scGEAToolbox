function [h] = i_scoreheatmap(Y, rowlabels, sce, parentfig)


    if nargin < 4, parentfig = []; end
    c = sce.c;
    
    hx=gui.myFigure;    
    h = heatmap(Y,'YDisplayLabels', rowlabels, ...
        'GridVisible','off','ColorbarVisible','off', ...
        'ColorScaling','scaledrows');
    h.XDisplayLabels = repmat({''}, 1, size(h.XDisplayData, 1));  % Remove X-axis labels

    %h.ColorScaling = 'scaledcolumns';  % Normalize each column
    %h.ColorScaling = 'scaledrows';    % Normalize each row
    %h.GridVisible = 'off';  % Turn off the borders around the cells
    %h.ColorbarVisible = 'off';
    % h = imagesc(Y);
    % set(gca, 'XTick', a-b);
    % set(gca, 'XTickLabel', cL);
    % set(gca, 'YTick', 1:length(glist));
    % set(gca, 'YTickLabel', glist);
    % set(gca, 'TickLength', [0, 0]);
    % box on
    %
    % szc = cumsum(szgn);
    % for k = 1:length(szc)
    %     xline(szc(k)+0.5, 'y-');
    % end

    hx.addCustomButton('on', {@gui.i_pickcolormap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
    hx.addCustomButton('off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
    hx.addCustomButton('on', @i_renamecat, 'guideicon.gif', 'Rename groups...');
    hx.addCustomButton('on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
    hx.addCustomButton('on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
    hx.addCustomButton('off', @i_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
    hx.addCustomButton('off', @i_flipxy, 'xplotpicker-geobubble2.gif', 'Flip XY');
    hx.show(parentfig);
    fliped = false;


    function i_flipxy(~, ~)
        %delete(h);
        fliped = ~fliped;
        if fliped
            h = imagesc(Y');
            set(gca, 'YTick', a-b);
            set(gca, 'YTickLabel', cL);
            %set(gca,'YTickLabelRotation',90);
            set(gca, 'XTick', 1:length(glist));
            set(gca, 'XTickLabel', glist);
            set(gca, 'XTickLabelRotation', 90);
            set(gca, 'TickLength', [0, 0]);
        else
            h = imagesc(Y);
            set(gca, 'XTick', a-b);
            set(gca, 'XTickLabel', cL);
            %set(gca,'XTickLabelRotation',0);
            set(gca, 'YTick', 1:length(glist));
            set(gca, 'YTickLabel', glist);
            set(gca, 'TickLength', [0, 0]);
        end
    end

    function i_renamecat(~, ~)
        tg = gui.i_inputgenelist(string(cL), true);
        if isempty(tg), return; end
        if length(tg) == length(cL)
            set(gca, 'XTick', a-b);
            set(gca, 'XTickLabel', tg(:))
            cL = tg;
        else
            errordlg('Wrong input.');
        end
    end

    function i_resetcolor(~, ~)
        set(gca, 'FontSize', 10);
        colormap default
    end

end
