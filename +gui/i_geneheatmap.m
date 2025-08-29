function [h] = i_geneheatmap(sce, thisc, glist, parentfig)

if nargin < 4, parentfig = []; end
if nargin < 3
    [glist] = gui.i_selectngenes(sce);
    if isempty(glist)
        gui.myHelpdlg(parentfig, 'No gene selected.');
        return;
    end
end
if nargin < 2
    [thisc, ~] = gui.i_select1class(sce,[],[],[],parentfig);
    if isempty(thisc), return; end
end
[c, cL, noanswer] = gui.i_reordergroups(thisc, [], parentfig);
if noanswer, return; end

[y, gidx] = ismember(upper(glist), upper(sce.g));
gidx = gidx(y);
glist = glist(y);

%[Xt]=gui.i_transformx(sce.X, [], [], parentfig);
%if isempty(Xt), return; end
Xt = sc_norm(sce.X);
Xt = log1p(Xt);

Y = Xt(gidx, :);
[~, cidx] = sort(c);
Y = Y(:, cidx);
[Y] = gui.i_norm4heatmap(Y);

% szgn = grpstats(c, c, @numel);
szgn = splitapply(@numel, c, c);
a = zeros(1, max(c));
b = zeros(1, max(c));
for k = 1:max(c)
    a(k) = sum(c <= k);
    b(k) = round(sum(c == k)./2);
end

% figure;
% heatmap(Y)
% assignin('base','Y',Y);
% assignin('base','g',glist) ;
% heatmap(Y,'YDisplayLabels',glist, ...
%     'XDisplayLabels',strings(size(Y,2),1), ...
%     'GridVisible',false,'ColorScaling','scaled',...
%     'ColorbarVisible',false)

hx=gui.myFigure(parentfig);
h = imagesc(Y);
% hFig.Colormap = repmat(linspace(0, 1, 25).', 1, 3);
set(gca, 'XTick', a-b);
set(gca, 'XTickLabel', cL);
%set(gca,'XTickLabelRotation',0);
set(gca, 'YTick', 1:length(glist));
set(gca, 'YTickLabel', glist);
set(gca, 'TickLength', [0, 0]);
% colormap(flipud(bone));
box on

szc = cumsum(szgn);
for k = 1:length(szc), xline(szc(k)+0.5, 'y-'); end

hx.addCustomButton('on', @in_callback_renamecat, 'guideicon.gif', 'Rename groups...');
hx.addCustomButton('off', @in_callback_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
hx.addCustomButton('off', @in_callback_flipxy, 'mat-wrap-text.jpg', 'Flip XY');

hx.show(parentfig);

fliped = false;

    function in_callback_flipxy(~, ~)
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

    function in_callback_renamecat(~, ~)
        tg = gui.i_inputgenelist(string(cL), true);
        if isempty(tg), return; end
        if length(tg) == length(cL)
            set(gca, 'XTick', a-b);
            set(gca, 'XTickLabel', tg(:))
            cL = tg;
        else
            gui.myErrordlg(parentfig, 'Wrong input.');
        end
    end

    function in_callback_resetcolor(~, ~)
        set(gca, 'FontSize', 10);
        colormap default
    end

end
