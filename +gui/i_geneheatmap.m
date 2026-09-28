function [h] = i_geneheatmap(sce, thisc, glist, parentfig)

if nargin < 4, parentfig = []; end
if nargin < 3
    [glist] = gui.i_selectngenes(sce, [], parentfig);
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

Xt = sc_norm(sce.X);
Xt = log1p(Xt);

Y = Xt(gidx, :);
[~, cidx] = sort(c);
Y = Y(:, cidx);
[Y] = gui.i_norm4heatmap(Y);

szgn = splitapply(@numel, c, c);
a = zeros(1, max(c));
b = zeros(1, max(c));
for k = 1:max(c)
    a(k) = sum(c <= k);
    b(k) = round(sum(c == k)./2);
end

hx=gui.myFigure(parentfig);
% Every draw below, including the button callbacks, targets this axes:
% after a dialog the current axes can be another window's.
ax = hx.AxHandle;
h = imagesc(ax, Y);
set(ax, 'XTick', a-b);
set(ax, 'XTickLabel', cL);
set(ax, 'YTick', 1:length(glist));
set(ax, 'YTickLabel', glist);
set(ax, 'TickLength', [0, 0]);
box(ax, 'on');

szc = cumsum(szgn);
for k = 1:length(szc), xline(ax, szc(k)+0.5, 'y-'); end

hx.addCustomButton('on', @in_callback_renamecat, 'guideicon.gif', 'Rename groups...');
hx.addCustomButton('off', @in_callback_resetcolor, 'plotpicker-geobubble2.gif', 'Reset color map');
hx.addCustomButton('off', @in_callback_flipxy, 'mat-wrap-text.jpg', 'Flip XY');

hx.show(parentfig);

fliped = false;

function in_callback_flipxy(~, ~)
        fliped = ~fliped;
        if fliped
            h = imagesc(ax, Y');
            set(ax, 'YTick', a-b);
            set(ax, 'YTickLabel', cL);
            set(ax, 'XTick', 1:length(glist));
            set(ax, 'XTickLabel', glist);
            set(ax, 'XTickLabelRotation', 90);
            set(ax, 'TickLength', [0, 0]);
        else
            h = imagesc(ax, Y);
            set(ax, 'XTick', a-b);
            set(ax, 'XTickLabel', cL);
            set(ax, 'YTick', 1:length(glist));
            set(ax, 'YTickLabel', glist);
            set(ax, 'TickLength', [0, 0]);
        end
    end

function in_callback_renamecat(~, ~)
        % The dialog belongs to this plot window, so closing it raises this
        % window rather than the main app.
        tg = gui.i_inputgenelist(string(cL), true, hx.FigHandle);
        if isempty(tg), return; end
        if length(tg) == length(cL)
            % Group labels sit on the Y axis while the map is flipped.
            if fliped
                set(ax, 'YTick', a-b, 'YTickLabel', tg(:));
            else
                set(ax, 'XTick', a-b, 'XTickLabel', tg(:));
            end
            cL = tg;
        else
            gui.myErrordlg(hx.FigHandle, 'Wrong input.');
        end
    end

function in_callback_resetcolor(~, ~)
        set(ax, 'FontSize', 10);
        colormap(ax, 'default')
    end

end
