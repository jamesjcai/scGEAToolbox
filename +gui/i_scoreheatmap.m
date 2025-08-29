function [h] = i_scoreheatmap(Y, rowlabels, sce, parentfig)


    if nargin < 4, parentfig = []; end
    
    [thisc, ~] = gui.i_select1class(sce,[],[],[],parentfig);
    if isempty(thisc), return; end
    [c, cL, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
    if noanswer, return; end


    % szgn = grpstats(c, c, @numel);
    szgn = splitapply(@numel, c, c);    
    a = zeros(1, max(c));
    b = zeros(1, max(c));
    for k = 1:max(c)
        a(k) = sum(c <= k);
        b(k) = round(sum(c == k)./2);
    end
    rowlabels = strrep(rowlabels, '_', '\_');
    hx=gui.myFigure(parentfig);    
    
    %{
    h = heatmap(Y,'YDisplayLabels', rowlabels, ...
        'GridVisible','off','ColorbarVisible','off', ...
        'ColorScaling','scaledrows');
    h.XDisplayLabels = repmat({''}, 1, size(h.XDisplayData, 1));  % Remove X-axis labels
    %h.ColorScaling = 'scaledcolumns';  % Normalize each column
    %h.ColorScaling = 'scaledrows';    % Normalize each row
    %h.GridVisible = 'off';  % Turn off the borders around the cells
    %h.ColorbarVisible = 'off';
    %}

    h = imagesc(Y);
    set(gca, 'XTick', a-b);
    set(gca, 'XTickLabel', cL);
    set(gca, 'YTick', 1:length(rowlabels));
    set(gca, 'YTickLabel', rowlabels);
    set(gca, 'TickLength', [0, 0]);
    box on

    szc = cumsum(szgn);
    for k = 1:length(szc)
        xline(szc(k)+0.5, 'y-');
    end

    hx.addCustomButton('off', @i_renamecat, 'edit.jpg', 'Rename groups...');
    hx.addCustomButton('off', @i_resetcolor, 'refresh.jpg', 'Reset color map');
    hx.addCustomButton('off', @i_flipxy, 'mat-wrap-text.jpg', 'Flip XY');
    hx.show(parentfig);
    fliped = false;

    function in_update
        if fliped
            h = imagesc(Y');
            set(gca, 'YTick', a-b);
            set(gca, 'YTickLabel', cL);
            %set(gca,'YTickLabelRotation',90);
            set(gca, 'XTick', 1:length(rowlabels));
            set(gca, 'XTickLabel', rowlabels);
            set(gca, 'XTickLabelRotation', 90);
            set(gca, 'TickLength', [0, 0]);
        else
            h = imagesc(Y);
            set(gca, 'XTick', a-b);
            set(gca, 'XTickLabel', cL);
            %set(gca,'XTickLabelRotation',0);
            set(gca, 'YTick', 1:length(rowlabels));
            set(gca, 'YTickLabel', rowlabels);
            set(gca, 'TickLength', [0, 0]);
        end
    end

    function i_flipxy(~, ~)
        %delete(h);
        fliped = ~fliped;
        in_update;
    end

    function i_renamecat(~, ~)
        tg = gui.i_inputgenelist(string(cL), true);
        if isempty(tg), return; end
        if length(tg) == length(cL)
            cL = tg;
            in_update;            
        else
            gui.myErrordlg(parentfig, 'Wrong input.');
        end
    end

    function i_resetcolor(~, ~)
        set(gca, 'FontSize', 10);
        colormap default
    end

end
