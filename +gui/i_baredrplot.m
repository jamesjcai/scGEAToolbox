function [hAx] = i_baredrplot(ax, c, t, parentfig, cellpicked)
% CELLPICKED, optional, is a logical mask over the points of AX's scatter;
% only those points, and the group labels on them, are kept in the copy,
% and the limits fit them.


hAx = [];
if nargin < 5, cellpicked = []; end
if nargin < 4, parentfig = []; end
if nargin < 3, t = 'tSNE'; end
if nargin < 2, c = []; end

colortag = 'k';

if ~isMATLABReleaseOlderThan('R2025a') && ~isempty(parentfig)
    try
        if strcmp('dark', parentfig.Theme.BaseColorStyle)
            colortag = [.5 .5 .5];
        end
    catch
        % keep default colortag if Theme.BaseColorStyle is unavailable
    end
end

isAxesHandle = isa(ax, 'matlab.graphics.axis.Axes'); % isgraphics(s, 'axes');
if ~isAxesHandle && isempty(c), error('Empty handle.'); end

hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;

if isAxesHandle
    if ~gui.i_isuifig(hFig)
        copyobj(ax.Children, hx.AxHandle);
        hAx = hx.AxHandle;
    else
        hAx = copyobj(ax, hFig);
    end
else
    if gui.i_isuifig(hFig)
        hAx = hx.AxHandle;


    else
        hAx = axes('Parent', hFig, 'Visible', 'off');
    end
end


set(hAx, 'XColor', 'none', 'YColor', 'none', 'ZColor', 'none');
grid(hAx, 'off');

% The group labels are data tips, and COPYOBJ brings them along on the
% copied scatter, so the copy needs no more. They used to be re-created on
% AX's scatter as well, which put a second set of labels on the main window
% each time this plot was opened.
if ~isempty(cellpicked) && ~all(cellpicked)
    in_keeppicked(cellpicked);
end

% Get axis limits and figure position
% ax = gca;
% ax = hAx;
xLimits = hAx.XLim;
yLimits = hAx.YLim;
zLimits = hAx.ZLim;

is3d1 = isprop(hAx, 'ZLim');
h = get(hAx, 'Children');
is3d2 = false;
if isscalar(h)
    if isprop(h, 'ZData')  % ismember('ZData', properties(h))
        is3d2 = any(arrayfun(@(x) ~isempty(get(x, 'ZData')), h));
    end
else
    for k=1:numel(h)
        if isa(h(k), 'matlab.graphics.chart.primitive.Scatter')
            is3d2 = any(arrayfun(@(x) ~isempty(get(x, 'ZData')), h(k)));
        end
    end
end
is3d = is3d1 & is3d2;

hold(hAx, 'on');

if is3d    % ======================================== 3D
    % Turn off the default axis display
    % set(gca, 'Visible', 'off');

    a = xLimits(1); b = yLimits(1); c = zLimits(1);
    la = xLimits(2)-a;
    lb = yLimits(2)-b;
    lc = zLimits(2)-c;

    % Draw custom arrows as axes using quiver3
    quiver3(hAx, a, b, c, la/5, 0, 0, 'Color', colortag, 'LineWidth', 1); % X-axis
    quiver3(hAx, a, b, c, 0, lb/5, 0, 'Color', colortag, 'LineWidth', 1); % Y-axis
    quiver3(hAx, a, b, c, 0, 0, lc/5, 'Color', colortag, 'LineWidth', 1); % Z-axis

    txt1 = sprintf('%s\\_1', t);
    txt2 = sprintf('%s\\_2', t);
    txt3 = sprintf('%s\\_3', t);

    % Label each arrow for clarity
    text(hAx, a+la/5, b, c, txt1);
    text(hAx, a, b+lb/5, c, txt2);
    text(hAx, a, b, c+lc/5, txt3);
    view(hAx, 3);
else          % ======================================== 2D
    hAx.Units = "pixels";
    r = hAx.Position(3)/hAx.Position(4);

    hAx.Units = "normalized";
    axPos = hAx.Position;

    % Convert data limits to figure normalized units
    % xArrowPos = [0, 1]; % normalized from left to right of the axes
    % yArrowPos = [0, 1]; % normalized from bottom to top of the axes

    % Draw x-axis arrow
    annotation(hFig, 'arrow', ...
        [axPos(1), axPos(1) + axPos(3)/7], ... % x positions
        [axPos(2), axPos(2)], ...            % y positions
        'Color', colortag, 'LineWidth', .5);

    % Draw y-axis arrow
    annotation(hFig, 'arrow', ...
        [axPos(1), axPos(1)], ...            % x positions
        [axPos(2), axPos(2) + r*(axPos(4)/7)], ... % y positions
        'Color', colortag, 'LineWidth', .5);

    txt1 = sprintf('%s\\_1', t);
    txt2 = sprintf('%s\\_2', t);

    textOpts = struct();
    textOpts.HorizontalAlignment = 'center';
    textOpts.VerticalAlignment = 'middle';
    textOpts.FontSize = 10;
    textOpts.FontWeight = 'normal';

    [~, b] = measureText(txt1, textOpts, hAx);
    text(hAx, xLimits(1), yLimits(1) - 2*b, txt1);

    [~, b] = measureText(txt2, textOpts, hAx);
    text(hAx, xLimits(1)-3*b, yLimits(1), txt2,'Rotation',90);
    view(hAx, 2);
end

hold(hAx,'off');
title(hAx, '')
subtitle(hAx, '')
xlim(hAx, xLimits);
ylim(hAx, yLimits);
hx.show(parentfig);


 function in_keeppicked(cellpicked)
    % Redraw the copied scatter with the picked points only. It is replaced
    % rather than trimmed: COPYOBJ leaves the copy sharing AX's scatter's
    % DataTipTemplate, so the copy's label rows cannot be changed without
    % changing the main window's, and the shared rows still describe every
    % cell.
    hold0 = ishold(hAx);
    hs = findobj(hAx, '-depth', 1, 'Type', 'scatter');
    if isempty(hs), return; end
    [~, imax] = max(arrayfun(@(x) numel(x.XData), hs));
    hs = hs(imax);
    npts = numel(hs.XData);
    assert(numel(cellpicked) == npts, ...
        'i_baredrplot: CELLPICKED must have one entry per plotted point.');
    cellpicked = cellpicked(:);

    % Labels on picked cells move to those cells' new positions.
    dts = findall(hs, 'Type', 'datatip');
    pos = zeros(numel(dts), 1);
    tiptext = cell(numel(dts), 1);
    for kd = 1:numel(dts)
        pos(kd) = dts(kd).DataIndex;
        tiptext{kd} = dts(kd).Content{1};
    end
    newindex = cumsum(cellpicked);
    keep = cellpicked(pos);
    pos = newindex(pos(keep));
    tiptext = tiptext(keep);

    % Per-point values are subset; a single colour or size stays as is.
    cdata = hs.CData;
    if size(cdata, 1) == npts, cdata = cdata(cellpicked, :); end
    sdata = hs.SizeData;
    if numel(sdata) == npts, sdata = sdata(cellpicked); end
    style = {'Marker', 'MarkerEdgeColor', 'MarkerFaceColor', ...
        'MarkerEdgeAlpha', 'MarkerFaceAlpha', 'LineWidth', 'Tag', ...
        'DisplayName', 'Clipping'};
    stylevalues = get(hs, style);

    hold(hAx, 'on');
    if isempty(hs.ZData)
        ns = scatter(hAx, hs.XData(cellpicked), hs.YData(cellpicked), ...
            sdata, cdata);
    else
        ns = scatter3(hAx, hs.XData(cellpicked), hs.YData(cellpicked), ...
            hs.ZData(cellpicked), sdata, cdata);
    end
    set(ns, style, stylevalues);
    delete(hs);
    if ~hold0, hold(hAx, 'off'); end

    stxtyes = cell(nnz(cellpicked), 1);
    stxtyes(pos) = tiptext;
    ns.DataTipTemplate.DataTipRows = dataTipTextRow('', stxtyes);
    for kd = 1:numel(pos)
        datatip(ns, 'DataIndex', pos(kd));
    end
    % Fit the picked cells rather than the whole embedding.
    set(hAx, 'XLimMode', 'auto', 'YLimMode', 'auto', 'ZLimMode', 'auto');
 end

 function [width, height] = measureText(txt, textOpts, ax)
    hTest = text(ax, 0, 0, txt, textOpts);
    textExt = get(hTest, 'Extent');
    delete(hTest);
    height = textExt(4)/3;    % Height
    width = textExt(3)/3;     % Width
 end

end
