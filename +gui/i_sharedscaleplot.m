function hFig = i_sharedscaleplot(values, glist, s, parentfig, viewAngle, valueLabel)
%I_SHAREDSCALEPLOT Draw genes on the embedding side by side, on one colour range.
%
%   hFig = gui.i_sharedscaleplot(values, glist, s, parentfig) draws one
%   panel per column of VALUES (cells x genes), named by GLIST, at the
%   embedding S, all on one colour range with one colour bar. Returns the
%   window.
%
%   gui.i_sharedscaleplot(..., viewAngle, valueLabel) sets the 3-D view
%   ([azimuth, elevation]; [] for the default) and the colour bar's label.
%
%   See also GUI.CALLBACK_SHOWGENESSHAREDSCALE.

arguments
    values {mustBeNumeric}
    glist string
    s {mustBeNumeric}
    parentfig = []
    viewAngle = []
    valueLabel = 'Expression'
end

numGenes = numel(glist);
is3d = size(s, 2) >= 3;
dotSize = 5;
cmapName = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');
isShared = true;
sharedRange = in_fullrange(values);

hx = gui.myFigure(parentfig);
hFig = hx.FigHandle;
hFig.Name = 'Genes on One Colour Scale';
hFig.Position(3) = max(hFig.Position(3), 380*numGenes);
scaleButton = hx.addCustomButton('off', @in_togglescale, 'setrange.gif', ...
    'Give each gene its own colour range');
hx.addCustomButton('off', @in_setrange, ...
    'gradient_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg', ...
    'Set the shared colour range...');

tl = tiledlayout(hFig, 1, numGenes, TileSpacing="compact", Padding="compact");
axs = gobjects(1, numGenes);
bars = gobjects(1, numGenes);
for k = 1:numGenes
    axs(k) = nexttile(tl);
    % High expression is drawn last, on top of the cells without it.
    [~, order] = sort(values(:, k), 'ascend');
    if is3d
        scatter3(axs(k), s(order, 1), s(order, 2), s(order, 3), dotSize, ...
            values(order, k), 'filled');
        if ~isempty(viewAngle)
            view(axs(k), viewAngle);
        end
    else
        scatter(axs(k), s(order, 1), s(order, 2), dotSize, values(order, k), ...
            'filled');
    end
    gui.i_setautumncolor(values(:), cmapName, true, true, axs(k), parentfig);
    title(axs(k), glist(k), Interpreter='none');
    subtitle(axs(k), gui.i_getsubtitle(values(:, k)));
    box(axs(k), 'on');
    bars(k) = colorbar(axs(k));
end
if is3d
    linkHandle = linkprop(axs, {'XLim', 'YLim', 'ZLim', 'View'});
else
    linkHandle = linkprop(axs, {'XLim', 'YLim'});
end
% LINKPROP stops linking once its object is cleared, so keep it with the
% window.
setappdata(hFig, 'SharedScaleLink', linkHandle);
if numGenes == 2
    detected = values > 0;
    title(tl, sprintf(['%s detected in %.1f%% of cells, %s in %.1f%%, ' ...
        'both in %.1f%%'], glist(1), 100*mean(detected(:, 1)), glist(2), ...
        100*mean(detected(:, 2)), 100*mean(all(detected, 2))), Interpreter='none');
end
in_applyscale();
hx.show(parentfig);

    function range = in_fullrange(v)
        % From zero, so that a colour reads as an absolute level; a
        % transform with negative values (residuals) extends it below.
        range = [min(0, min(v(:))), max(v(:))];
        if ~(range(2) > range(1))
            range = range(1) + [0, 1];
        end
    end

    function in_applyscale()
        for kk = 1:numGenes
            if isShared
                clim(axs(kk), sharedRange);
                % One bar, on the last panel, serves every panel.
                bars(kk).Visible = kk == numGenes;
            else
                clim(axs(kk), in_fullrange(values(:, kk)));
                bars(kk).Visible = true;
            end
            bars(kk).Label.String = '';
        end
        if isShared
            bars(numGenes).Label.String = sprintf('%s, shared range', valueLabel);
        else
            bars(numGenes).Label.String = valueLabel;
        end
    end

    function in_togglescale(~, ~)
        isShared = ~isShared;
        if isShared
            scaleButton.TooltipString = 'Give each gene its own colour range';
        else
            scaleButton.TooltipString = 'Put every gene on one colour range';
        end
        in_applyscale();
    end

    function in_setrange(~, ~)
        answer = gui.myInputdlg({'Lowest value:', 'Highest value:'}, ...
            'Shared Colour Range', {sprintf('%.4g', sharedRange(1)), ...
            sprintf('%.4g', sharedRange(2))}, hFig);
        if isempty(answer)
            return;
        end
        range = str2double(string(answer(:))).';
        if any(isnan(range)) || ~(range(2) > range(1))
            gui.myWarndlg(hFig, ['Enter two numbers, the highest larger ' ...
                'than the lowest.']);
            return;
        end
        sharedRange = range;
        isShared = true;
        scaleButton.TooltipString = 'Give each gene its own colour range';
        in_applyscale();
    end
end
