function [hFig] = i_dotplot(X, g, c, cL, tgene, uselog, ...
    ttxt, parentfig, varargin)
% I_DOTPLOT  Seurat-style dot plot: dot size = fraction of cells in a group
% expressing the gene, dot colour = average expression in that group.
%
% Name-value options (after PARENTFIG):
%   'InputType' - what X holds, which decides how the average is taken:
%       'auto'   (default) non-negative integers -> 'counts'; other
%                non-negative values -> 'log'; any negative -> 'other'.
%       'counts' library-size normalise to 1e4 per cell (library size from
%                all rows of X, so pass the full matrix), then
%                log1p(mean(x)).
%       'log'    X is log1p of normalised data: log1p(mean(expm1(x))), i.e.
%                average in linear space, as Seurat's DotPlot does.
%       'other'  plain mean(x); USELOG then applies log1p to it (the
%                behaviour of this function before 'InputType' existed).
%   'Scale'     - true (default) z-scores each gene across groups and clips
%                to [-2.5 2.5], so colour compares groups within a gene;
%                false keeps the unscaled averages on one shared scale.

p = inputParser;
addParameter(p, 'InputType', 'auto', @(s) any(strcmpi(s, ...
    {'auto', 'counts', 'log', 'other'})));
addParameter(p, 'Scale', true, @(v) islogical(v) || isnumeric(v));
parse(p, varargin{:});
inputtype = lower(p.Results.InputType);
doscale = logical(p.Results.Scale);

if nargin < 8, parentfig = []; end
if nargin < 7, ttxt = []; end
DOTSIZE = 0.5;
cL=cL(:);

if nargin < 6 || isempty(uselog), uselog = true; end
[yes, gidx] = ismember(tgene, g);
if ~any(yes)
    warning('No genes found.');
    return;
end
z = length(tgene) - sum(yes);
if z > 0
    fprintf('%d gene(s) not in the list are excluded.\n', z);
end
tgene = tgene(yes);
gidx = gidx(yes);

Xs = full(X(gidx, :));
if strcmp(inputtype, 'auto')
    if any(Xs(:) < 0)
        inputtype = 'other';
    elseif all(mod(Xs(:), 1) == 0)
        inputtype = 'counts';
    else
        inputtype = 'log';
    end
end
switch inputtype
    case 'counts'
        libsize = full(sum(X, 1));
        libsize(libsize == 0) = 1;
        Xs = Xs ./ libsize * 1e4;
    case 'log'
        Xs = expm1(Xs);
end

nG = length(tgene);
nC = length(cL);
D = zeros(nG, nC);          % average expression (log scale unless 'other')
P = zeros(nG, nC);          % fraction of cells expressing
for kk = 1:nC
    a0 = Xs(:, c == kk);
    D(:, kk) = mean(a0, 2);
    P(:, kk) = mean(a0 ~= 0, 2);
end
if ~strcmp(inputtype, 'other') || uselog
    D = log1p(D);
end

if doscale
    Dz = (D - mean(D, 2)) ./ std(D, 0, 2);
    Dz(~isfinite(Dz)) = 0;  % a gene flat across groups, or a single group
    Dz = min(max(Dz, -2.5), 2.5);
else
    Dz = D;
end

% Long format, groups varying fastest within each gene
x = repmat((1:nC)', nG, 1);
y = repelem((1:nG)', nC);
AvgExpr = reshape(D.', [], 1);
AvgExprScaled = reshape(Dz.', [], 1);
PrtExpr = reshape(P.', [], 1);
sz = PrtExpr;
vl = AvgExprScaled;

GroupList = repmat(string(cL), length(tgene), 1);
GeneList = [];
for kx = 1:length(tgene)
    GeneList = [GeneList; repmat(tgene(kx), length(cL), 1)];
end

T = table(GeneList, GroupList, AvgExpr, AvgExprScaled, PrtExpr);


txgene = [" "; tgene(:)];

hx=gui.myFigure(parentfig);
hFig = hx.FigHandle;
% Every draw below, the button callbacks included, targets this axes: the
% current axes after a dialog can be another window's.
ax = hx.AxHandle;

dotsz = DOTSIZE;
sz(sz == 0) = eps;
vl = vl + 0.001;
afa = scatter(ax, x, y, dotsz*500*sz, vl, 'filled');
hold(ax, 'on');
afb = scatter(ax, x, y, dotsz*500*sz, 'k');

af{1} = scatter(ax, max(x)+1, 1, dotsz*500*0.75, 'k');
af{2} = text(ax, max(x)+1.4, 1, '75%', 'BackgroundColor', 'none');
af{3} = scatter(ax, max(x)+1, 2, dotsz*500*0.5, 'k');
af{4} = text(ax, max(x)+1.4, 2, '50%', 'BackgroundColor', 'none');
af{5} = scatter(ax, max(x)+1, 3, dotsz*500*0.25, 'k');
af{6} = text(ax, max(x)+1.4, 3, '25%', 'BackgroundColor', 'none');


xmax0 = length(cL) + 2.5;   % default right x limit; widened to fit legend labels
xlim(ax, [0.5, xmax0]);
ylim(ax, [0.5, max([4, length(txgene)]) - 0.5]);
% colorbar
% colorbar('northoutside');

set(ax, 'YTick', 0:length(tgene))
set(ax, 'YTickLabel', txgene)
set(ax, 'XTick', 0:length(cL))
cL = gui.i_escapeunderscore(cL);
set(ax, 'XTickLabel', [{''}; cL(:); {''}])
colormap(ax, flipud(summer));
box(ax, 'on');
grid(ax, 'on');

cb = colorbar(ax, 'eastoutside');
if doscale
    cb.Label.String = 'Avg. expr. (scaled)';
else
    cb.Label.String = 'Avg. expr.';
end
axposition = ax.Position;
cb.Position(3) = cb.Position(3) * 0.5;
cb.Position(4) = cb.Position(4) * min(1, 5 / length(tgene));
ax.Position = axposition;
in_fitlegend(getpixelposition(ax) * [0; 0; 1; 0]);

if ~isempty(ttxt)
    ttxt = gui.i_escapeunderscore(ttxt);
    title(ttxt);
end
hx.addCustomButton('off', @i_renamecat, 'edit.jpg', 'Rename groups...');
hx.addCustomButton('off', @in_callback_savetable, 'floppy-disk-arrow-in.jpg', 'Export data...');
hx.addCustomButton('off', @in_callback_resetcolor, 'refresh_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', 'Reset Colormap');
sqTip = 'Square Grid (resize window for equal x/y spacing)';
sqSaved = [];               % layout to restore when the toggle is released
hx.addCustomButton('off', @in_callback_squaregrid, 'view-grid.jpg', sqTip);
% myFigure only makes push tools; swap this one for a toggle tool in the
% same place so its pressed state shows whether the square grid is on.
pt = findall(hFig, 'Type', 'uipushtool', 'Tooltip', sqTip);
sqBtn = uitoggletool(pt.Parent, 'CData', pt.CData, 'Tooltip', sqTip, ...
    'Separator', pt.Separator, 'ClickedCallback', @in_callback_squaregrid);
delete(pt);
% gui.myFigure's Enlarge Window button and the square grid both resize the
% window, so each releases the other first: otherwise the second one to be
% undone restores a layout that is stale. The button stays where myFigure
% put it; only what a click does is taken over here.
winBtn = hx.SizeButton;
winBtn.ClickedCallback = @in_callback_winsize;
if nargout > 0, return; end
hx.show(parentfig);


function in_callback_savetable(~, ~)
        answer = gui.myQuestdlg(hFig, 'Export & save data to:', '', ...
            {'Workspace', 'TXT/CSV file', 'Excel file'}, 'Workspace');
        if ~isempty(answer)
            GroupList = repmat(string(cL), length(tgene), 1);
            GeneList = [];
            for k = 1:length(tgene)
                GeneList = [GeneList; repmat(tgene(k), length(cL), 1)];
            end
            T = table(GeneList, GroupList, AvgExpr, AvgExprScaled, PrtExpr);
            switch answer
                case 'Workspace'
                    labels = {'Save T to variable named:', 'Save D to variable named:'};
                    vars = {'T', 'D'};
                    values = {T, D};
                    [~, ~] = export2wsdlg(labels, vars, values, ...
                        'Save Data to Workspace');
                case 'TXT/CSV file'
                    [file, path] = uiputfile({'*.csv'; '*.*'}, 'Save as');
                    % Back to this plot window, not the main app.
                    gui.i_raisefig(hFig);
                    if isequal(file, 0) || isequal(path, 0)
                        return;
                    else
                        fw = gui.myWaitbar(hFig);
                        filename = fullfile(path, file);
                        writetable(T, filename, 'FileType', 'text');
                        OKPressed = true;
                        gui.myWaitbar(hFig, fw);
                    end
                case 'Excel file'

                    [file, path] = uiputfile({'*.xlsx'; '*.*'}, 'Save as');
                    % Back to this plot window, not the main app.
                    gui.i_raisefig(hFig);
                    if isequal(file, 0) || isequal(path, 0)
                        return;
                    else
                        fw = gui.myWaitbar(hFig);
                        filename = fullfile(path, file);
                        writetable(T, filename, 'FileType', 'spreadsheet');
                        OKPressed = true;
                        gui.myWaitbar(hFig, fw);
                    end
            end
        end
    end

function in_callback_resizedot(~, ~)
        dotsz = dotsz * 0.9;
        if dotsz < 0.2, dotsz = 1.0; end
        delete(afa);
        delete(afb);
        delete(af{1});
        delete(af{3});
        delete(af{5});
        afa = scatter(ax, x, y, dotsz*500*sz, vl, 'filled');
        hold(ax, 'on');
        afb = scatter(ax, x, y, dotsz*500*sz, 'k');
        af{1} = scatter(ax, max(x)+1, 1, dotsz*500*0.75, 'k');
        af{3} = scatter(ax, max(x)+1, 2, dotsz*500*0.5, 'k');
        af{5} = scatter(ax, max(x)+1, 3, dotsz*500*0.25, 'k');
    end

function i_renamecat(~, ~)
        tg = gui.i_inputgenelist(string(cL), true, hFig);
        if isempty(tg), return; end
        if length(tg) == length(cL)
            set(ax, 'XTick', 0:length(cL));
            set(ax, 'XTickLabel', [{''}; tg(:); {''}])
            cL = tg;
        else
            gui.myErrordlg(hFig, 'Wrong input.');
        end
    end

function in_callback_squaregrid(src, ~)
        % Toggle: pressed squares the dot grid, released puts the window and
        % axes back the way they were before it was pressed.
        if strcmp(src.State, 'on')
            if strcmp(hFig.WindowStyle, 'docked')
                gui.myHelpdlg(hFig, 'Undock the figure to resize it.', '');
                src.State = 'off';
                return;
            end
            if ~isempty(winBtn.UserData)            % enlarged: restore first
                gui.i_togglewinsize(winBtn, [], hFig);
                in_setfigpos(hFig.Position);
                in_fitlegend(getpixelposition(afa.Parent) * [0; 0; 1; 0]);
            end
            sqSaved = in_getlayout();
            [fpos, axT, cbT, xlim2] = in_squarelayout();
            in_applylayout(fpos, axT, cbT, xlim2);
            src.Tooltip = 'Restore Window Size';
        else
            if ~isempty(sqSaved)
                in_applylayout(sqSaved{:});
                sqSaved = [];
            end
            src.Tooltip = sqTip;
        end
    end

function in_callback_winsize(src, event)
        if strcmp(sqBtn.State, 'on')                % square grid: release it
            sqBtn.State = 'off';
            in_callback_squaregrid(sqBtn);
        end
        gui.i_togglewinsize(src, event, hFig);
        % The size-legend labels have a fixed pixel width, so the right x
        % limit that fits them depends on how wide the axes now are.
        in_setfigpos(hFig.Position);
        in_fitlegend(getpixelposition(afa.Parent) * [0; 0; 1; 0]);
    end

function L = in_getlayout()
        % Window, axes and colorbar in pixels, plus the right x limit.
        ax1 = afa.Parent;
        cb1 = findobj(hFig, 'Type', 'colorbar');
        cbT = [];
        if ~isempty(cb1), cbT = getpixelposition(cb1(1)); end
        L = {getpixelposition(hFig), getpixelposition(ax1), cbT, ax1.XLim(2)};
    end

function [fpos, axT, cbT, xlim2] = in_squarelayout()
        % Layout in which one x unit spans as many pixels as one y unit, so
        % the dots sit on a square grid. The coarser axis spacing is kept
        % and the other axis grows to match, capped so the window fits on
        % screen. Margins (tick labels, title, colorbar) keep their pixel
        % size. The right x limit leaves room for the legend labels, a fixed
        % pixel width, so it depends on the spacing and is solved with it.
        ax1 = afa.Parent;
        cb1 = findobj(hFig, 'Type', 'colorbar');
        drawnow;
        fpos = getpixelposition(hFig);
        pa = getpixelposition(ax1);
        mL = pa(1); mB = pa(2);
        mR = fpos(3) - pa(1) - pa(3);
        mT = fpos(4) - pa(2) - pa(4);
        x1 = ax1.XLim(1);
        yr = diff(ax1.YLim);
        [tw, xt] = in_legendtext(ax1);
        xr = @(u) max(xmax0, xt + (tw + 6) / u) - x1;
        scr = get(groot, 'ScreenSize');
        availW = 0.95*scr(3) - mL - mR;
        availH = 0.90*scr(4) - mB - mT;
        u = max(pa(3) / diff(ax1.XLim), pa(4) / yr);     % pixels per unit
        u = min(u, availH / yr);
        for iter = 1:20                                   % u*xr(u) grows with u
            if u * xr(u) <= availW, break; end
            u = u * 0.95 * availW / (u * xr(u));
        end
        axW = u * xr(u);
        axH = u * yr;
        newW = mL + axW + mR;
        newH = mB + axH + mT;
        top = fpos(2) + fpos(4);                          % keep the title bar in place
        fpos = [min(max(fpos(1), 1), scr(3) - newW), max(top - newH, 1), newW, newH];
        axT = [mL, mB, axW, axH];
        xlim2 = x1 + xr(u);
        cbT = [];
        if ~isempty(cb1)
            pc = getpixelposition(cb1(1));
            cbT = [mL + axW + pc(1) - (pa(1) + pa(3)), mB + pc(2) - pa(2), pc(3), pc(4)];
        end
    end

function in_applylayout(fpos, axT, cbT, xlim2)
        % Resize the window and place axes and colorbar at pixel targets.
        ax1 = afa.Parent;
        cb1 = findobj(hFig, 'Type', 'colorbar');
        if isempty(cbT), cb1 = []; end
        restoreunits = in_pixelunits(ax1, cb1);
        in_setfigpos(fpos);
        ax1.Position = axT;
        ax1.XLim(2) = xlim2;
        if ~isempty(cb1), cb1(1).Position = cbT; end
        restoreunits();
        % hFig.Position reports the requested size before the OS applies
        % it, so the units restore above can convert against the old size
        % (seen once after exportgraphics). Re-apply pixel targets until
        % they hold.
        for iter = 1:20
            drawnow;
            drift = any(abs(getpixelposition(ax1) - axT) > 1);
            if ~isempty(cb1)
                drift = drift || any(abs(getpixelposition(cb1(1)) - cbT) > 1);
            end
            if ~drift, break; end
            setpixelposition(ax1, axT);
            if ~isempty(cb1), setpixelposition(cb1(1), cbT); end
            pause(0.05);
        end
    end

function in_fitlegend(axW)
        % Widen the right x limit just enough that the size-legend labels,
        % which have a fixed pixel width, end inside the axes box for an axes
        % axW pixels wide. Solves (x2 - xt) * axW / (x2 - x1) = tw + 6.
        ax1 = afa.Parent;
        [tw, xt] = in_legendtext(ax1);
        Tpx = tw + 6;
        x1 = ax1.XLim(1);
        if axW <= Tpx, return; end
        ax1.XLim(2) = max(xmax0, (xt * axW - Tpx * x1) / (axW - Tpx));
    end

function [tw, xt] = in_legendtext(ax1)
        % Widest size-legend label in pixels, and its left edge in data x.
        tw = 0;
        xt = ax1.XLim(2);
        for k = [2 4 6]
            if ~isvalid(af{k}), continue; end
            af{k}.Units = 'pixels';
            tw = max(tw, af{k}.Extent(3));
            af{k}.Units = 'data';
            xt = min(xt, af{k}.Position(1));
        end
    end

function restorefn = in_pixelunits(ax1, cb1)
        % Pin figure, axes and colorbar to pixel units so a window resize
        % leaves their size and offsets alone; returns a restore function.
        drawnow;
        h = [hFig; ax1; cb1(:)];
        old = get(h, 'Units');
        set(h, 'Units', 'pixels');
        restorefn = @() set(h, {'Units'}, cellstr(old));
    end

function in_setfigpos(fpos)
        % The OS applies a resize asynchronously; wait for it to land before
        % anything is placed against the new size.
        hFig.Position = fpos;
        for iter = 1:20
            drawnow;
            if all(abs(hFig.Position(3:4) - fpos(3:4)) < 2), break; end
            pause(0.05);
        end
    end

function in_callback_resetcolor(~, ~)
        dotsz = DOTSIZE;
        set(ax, 'FontSize', 10);
        in_callback_resizedot;
        colormap(ax, flipud(summer));
    end

end
