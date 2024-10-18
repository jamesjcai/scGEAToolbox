function [h] = i_scoreheatmap(Y, rowlabels, sce, parentfig)


if nargin < 4, parentfig = []; end
c = sce.c;
% if nargin < 2
%     [thisc, ~] = gui.i_select1class(sce);
%     if isempty(thisc), return; end
% end
% [c, cL, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
% if noanswer, return; end
% [~, cidx] = sort(c);
% Y = Y(:, cidx);
% [Y] = gui.i_norm4heatmap(Y);

% szgn = grpstats(c, c, @numel);
% a = zeros(1, max(c));
% b = zeros(1, max(c));
% for k = 1:max(c)
%     a(k) = sum(c <= k);
%     b(k) = round(sum(c == k)./2);
% end

% figure;
% heatmap(Y)
% assignin('base','Y',Y);
% assignin('base','g',glist) ;
% heatmap(Y,'YDisplayLabels',glist, ...
%     'XDisplayLabels',strings(size(Y,2),1), ...
%     'GridVisible',false,'ColorScaling','scaled',...
%     'ColorbarVisible',false)

f = figure('Visible', 'off');

if ~isempty(parentfig)
    p = parentfig.Position;
    cx = [p(1)+p(3)/2 p(2)+p(4)/2];
    px = f.Position;
    px_new = [cx(1)-px(3)/2 cx(2)-px(4)/2];
    movegui(f, px_new);
else
    movegui(f, 'center');
end


    
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

% szc = cumsum(szgn);
% for k = 1:length(szc)
%     xline(szc(k)+0.5, 'y-');
% end
tb = uitoolbar('Parent', f);
pkg.i_addbutton2fig(tb, 'on', {@gui.i_pickcolormap, c}, 'plotpicker-compass.gif', 'Pick new color map...');
pkg.i_addbutton2fig(tb, 'off', @gui.i_changefontsize, 'noun_font_size_591141.gif', 'ChangeFontSize');
pkg.i_addbutton2fig(tb, 'on', @i_renamecat, 'guideicon.gif', 'Rename groups...');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
pkg.i_addbutton2fig(tb, 'off', @i_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
pkg.i_addbutton2fig(tb, 'off', @i_flipxy, 'xplotpicker-geobubble2.gif', 'Flip XY');

set(f, 'visible', 'on');
drawnow;

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
