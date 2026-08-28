function [h] = i_scoreheatmap(Y, rowlabels, sce, parentfig, thisc)
% THISC, when given, is the grouping already chosen by the caller - one
% label per column of Y. Passing it in matters for more than saving a
% click: the caller may have grouped by several crossed variables, which
% the i_select1class dialog below has no way to offer, so re-asking there
% would quietly describe a different partition than the one the scores
% were compared under.

if nargin < 4, parentfig = []; end
if nargin < 5 || isempty(thisc)
    [thisc, ~] = gui.i_select1class(sce,[],[],[],parentfig);
end
if isempty(thisc), return; end

% Skip the reorder prompt when there are too many groups for it to be
% usable - i_reordergroups insists every item be reselected by hand, which
% is a hostile ask at dozens of levels. Fall back to its natural order.
if numel(unique(thisc)) > 30
    [c, cL] = findgroups(string(thisc));
else
    [c, cL, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
    if noanswer, return; end
end

% The tick positions and divider lines below are cumulative group sizes, so
% they only describe the image if its columns are actually in group order.
% Y arrives in original cell order, so sort both together - without this
% the labels and the yellow separators mark boundaries the picture does not
% have.
[c, ord] = sort(c(:));
Y = Y(:, ord);

% szgn = grpstats(c, c, @numel);
szgn = splitapply(@numel, c, c);
a = zeros(1, max(c));
b = zeros(1, max(c));
for k = 1:max(c)
    a(k) = sum(c <= k);
    b(k) = round(sum(c == k)./2);
end
rowlabels = gui.i_escapeunderscore(rowlabels);
hx=gui.myFigure(parentfig);

%{
h = heatmap(Y,'YDisplayLabels', rowlabels, ...
    'GridVisible','off','ColorbarVisible','off', ...
    'ColorScaling','scaledrows');
h.XDisplayLabels = repmat({''}, 1, size(h.XDisplayData, 1));  % Remove X-axis labels
% h.ColorScaling = 'scaledcolumns';  % Normalize each column
% h.ColorScaling = 'scaledrows';    % Normalize each row
% h.GridVisible = 'off';  % Turn off the borders around the cells
% h.ColorbarVisible = 'off';
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
    % set(gca,'YTickLabelRotation',90);
    set(gca, 'XTick', 1:length(rowlabels));
    set(gca, 'XTickLabel', rowlabels);
    set(gca, 'XTickLabelRotation', 90);
    set(gca, 'TickLength', [0, 0]);
else
    h = imagesc(Y);
    set(gca, 'XTick', a-b);
    set(gca, 'XTickLabel', cL);
    % set(gca,'XTickLabelRotation',0);
    set(gca, 'YTick', 1:length(rowlabels));
    set(gca, 'YTickLabel', rowlabels);
    set(gca, 'TickLength', [0, 0]);
end
end

function i_flipxy(~, ~)
% delete(h);
fliped = ~fliped;
in_update;
end

function i_renamecat(~, ~)
tg = gui.i_inputgenelist(string(cL), true, parentfig);
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
