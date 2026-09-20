function callback_CrossTabulation(src, ~)

[FigureHandle, sce] = gui.gui_getfigsce(src);

% [thisc1, clabel1, thisc2, clabel2] = gui.i_select2class(sce, true, FigureHandle);
[thisc1, clabel1, thisc2, clabel2] = gui.i_select2states(sce, true, FigureHandle);

% [answer] = gui.myQuestdlg(FigureHandle, 'Manually order groups?', '');
if isempty(thisc1), return; end

iscon1 = i_iscontinuous(thisc1);
iscon2 = ~isempty(thisc2) && i_iscontinuous(thisc2);

if ~isempty(thisc2) && (iscon1 || iscon2)
    if iscon1 && iscon2
        % Both continuous: scatter plot
        hx = gui.myFigure(FigureHandle);
        axbar = hx.AxHandle;
        if isempty(axbar), axbar = gca; end
        scatter(ax, thisc1, thisc2, 5, 'filled', 'MarkerFaceAlpha', 0.3);
        xlabel(axbar, strrep(clabel1, '_', '\_'));
        ylabel(axbar, strrep(clabel2, '_', '\_'));
        title(ax, sprintf('%s vs %s', clabel1, clabel2));
        [rho, pval] = corr(double(thisc1(:)), double(thisc2(:)), ...
            'Type', 'Spearman', 'Rows', 'complete');
        subtitle(axbar, sprintf('Spearman \\rho=%.3f, p=%.2e', rho, pval));
        hx.show(FigureHandle);
        return;
    end
    % One continuous, one categorical
    if iscon1
        conc = double(thisc1(:));
        conlabel = clabel1;
        catc = thisc2;
        catlabel = clabel2;
    else
        conc = double(thisc2(:));
        conlabel = clabel2;
        catc = thisc1;
        catlabel = clabel1;
    end
    gui.i_violinplot(conc, catc, ...
        sprintf('%s by %s', conlabel, catlabel), ...
        true, [], [], FigureHandle);
    return;
end

[c, cL1, noanswer] = gui.i_reordergroups(thisc1, [], FigureHandle);
if noanswer, return; end
if isnumeric(thisc1)
    thisc1 = categorical(cL1(c), cL1);
else
    thisc1 = categorical(thisc1, cL1);
end
if ~isempty(thisc2)
    [c2, cL2, noanswer] = gui.i_reordergroups(thisc2, [], FigureHandle);
    if noanswer, return; end
    if isnumeric(thisc2)
        thisc2 = categorical(cL2(c2), cL2);
    else
        thisc2 = categorical(thisc2, cL2);
    end
else
    cL2 = [];
end

if isempty(thisc1), return; end
if isempty(thisc2) || isempty(cL2)

    hx = gui.myFigure(FigureHandle);
    Tb = tabulate(thisc1);
    y0 = Tb(:,2);
    if iscell(y0), y0 = cell2mat(y0); end
    y0 = y0(:);
    labels0 = string(Tb(:,1));
    axbar = hx.AxHandle;
    if isempty(axbar), axbar = gca; end

    % 0 = the order the groups came in, which GUI.I_REORDERGROUPS may have
    % had the user set by hand, so it has to stay reachable rather than be
    % overwritten by the first sort.
    sortstate = 0;
    ispie = false;
    in_drawonebar();
    hx.addCustomButton('off', @in_callback_sortbars, "reorder.jpg", ...
        'Sort groups by cell count');
    hx.addCustomButton('off', @in_callback_togglepie, "plotpicker-pie.gif", ...
        'Switch between bar and pie chart');
    hx.show(FigureHandle);
    return;
end


fw = gui.myWaitbar(FigureHandle);

% sizesorted = false;
labelsx='';
labelsy='';
T=[];

hx = gui.myFigure(FigureHandle, true);
tab=cell(2,1);
ax0=cell(2,1);
ax=cell(2,2);

% Each tab crosstabs the two variables the other way round, so it has its
% own table, its own x groups and its own legend. Held per tab rather than
% in one shared set of variables, so redrawing one tab in a new order
% cannot plot it with the other tab's labels.
Tall = cell(2,1);
labelsxall = cell(2,1);
labelsyall = cell(2,1);
clabelall = cell(2,1);
llabelall = cell(2,1);

% 0 = the order the groups came in, which GUI.I_REORDERGROUPS may have had
% the user set by hand. 1 = descending by cell count, 2 = ascending.
sortstate = 0;
ismosaic = false;

tabgp = uitabgroup();
for k=1:2
    switch k
        case 1
            thiscA = thisc1;
            thiscB = thisc2;
            clabel = clabel1;
            llabel = clabel2;
        case 2
            thiscA = thisc2;
            thiscB = thisc1;
            clabel = clabel2;
            llabel = clabel1;
    end
    in_crossplot(thiscA, thiscB);
    Tall{k} = T;
    labelsxall{k} = labelsx;
    labelsyall{k} = labelsy;
    clabelall{k} = clabel;
    llabelall{k} = llabel;
    tab{k} = uitab(tabgp, 'Title', sprintf('Tab%d',k));
    % tab{k} = uitab(tabgp, 'Title', sprintf('%s-%s',clabel,llabel));
    ax0{k} = axes('parent',tab{k});
    ax{k,1} = subplot(2,1,1);
    ax{k,2} = subplot(2,1,2);
    in_drawtab(k);
end

hx.addCustomButton('off', @in_callback_saveCrossTable, "floppy-disk-arrow-in.jpg", 'Save cross-table');
hx.addCustomButton('off', @in_callback_sortgroups, "reorder.jpg", 'Sort groups by cell count');
hx.addCustomButton('off', @in_callback_togglemosaic, ...
    "full_stacked_bar_chart_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", ...
    'Mosaic: bar width by group size');
gui.myWaitbar(FigureHandle, fw);
hx.show(FigureHandle);


function in_drawonebar()
        % Redraw the one-variable plot in the current sort order, as bars or
        % as a pie. Both read the same counts, so the sort carries across.
        switch sortstate
            case 1
                [ybar, idxbar] = sort(y0, 'descend');
                statetxt = 'sorted by count, descending';
            case 2
                [ybar, idxbar] = sort(y0, 'ascend');
                statetxt = 'sorted by count, ascending';
            otherwise
                ybar = y0;
                idxbar = (1:numel(y0))';
                statetxt = '';
        end
        labelsbar = labels0(idxbar);

        % A pie sets axis equal and hides the ruler, and bars need the ruler
        % back, so reset rather than clear. The legend is not a child of the
        % axes and would outlive either.
        legend(axbar, 'off');
        cla(axbar, 'reset');
        colormap(axbar, "turbo")
        if ispie
            in_drawpie(ybar, labelsbar);
        else
            in_drawbars(ybar, labelsbar);
        end
        subtitle(axbar, statetxt);
    end

function in_drawbars(yvals, labs)
        % Counts and labels come in as parameters: they live only inside
        % these sibling nested functions, never in the enclosing one, so
        % MATLAB gives each its own copy rather than sharing them.
        bar(axbar, yvals, 'FaceColor', "flat");
        xticks(axbar, 1:numel(labs));
        xticklabels(axbar, strrep(labs, '_', '\_'));
        ylabel(axbar, '# of cells');
        xlabel(axbar, strrep(clabel1, '_', '\_'));
    end

function in_drawpie(yvals, labs)
        % PIE drops zero and negative slices, which would slide the labels
        % onto the wrong wedges, so drop them here instead.
        keep = yvals > 0;
        ndropped = sum(~keep);
        ypie = yvals(keep);
        labelspie = labs(keep);
        if isempty(ypie)
            title(axbar, 'No cells to show');
            return;
        end

        pie(axbar, ypie);
        % The slice names go in the legend rather than on the wedges: group
        % names here are cell types, too long to sit beside a thin slice.
        % Counts ride along, since the wedges themselves only show percent.
        keytxt = strrep(labelspie, '_', '\_') + " (" + string(ypie) + ")";
        legend(axbar, keytxt, 'Location', 'eastoutside');
        ttl = strrep(clabel1, '_', '\_');
        if ndropped > 0
            ttl = ttl + sprintf(' (%d empty group(s) not shown)', ndropped);
        end
        title(axbar, ttl);
    end

function in_callback_sortbars(~, ~)
        % Descending, then ascending, then back to the order it started in.
        sortstate = mod(sortstate+1, 3);
        in_drawonebar();
    end

function in_callback_togglepie(~, ~)
        ispie = ~ispie;
        in_drawonebar();
    end

function in_crossplot(thiscA, thiscB)
        % t = table(thiscA, thiscB);
        % t = sortrows(t, [1, 2]);
        % thiscA = t.thiscA;
        % thiscB = t.thiscB;

        % iscategorical(thiscA)
        % iscategorical(thiscB)
        % categories(thiscA)
        % categories(thiscB)

        thiscB = reordercats(thiscB, flipud(categories(thiscB)));

        [T, ~, ~, labelsxy] = crosstab(thiscA, thiscB);

        labelsx = labelsxy(:, 1);
        labelsx = labelsx(~cellfun('isempty', labelsx));
        labelsy = labelsxy(:, 2);
        labelsy = labelsy(~cellfun('isempty', labelsy));

    end


function in_drawtab(kdraw)
        % Redraw both panels of one tab in the current sort order.
        [idx, statetxt] = in_grouporder(kdraw);
        in_plot1(kdraw, idx, statetxt);
        in_plot2(kdraw, idx);
    end

function [idx, statetxt] = in_grouporder(korder)
        % Row order for one tab: its groups by how many cells they hold.
        rowtotal = sum(Tall{korder}, 2);
        switch sortstate
            case 1
                [~, idx] = sort(rowtotal, 'descend');
                statetxt = 'groups sorted by cell count, descending';
            case 2
                [~, idx] = sort(rowtotal, 'ascend');
                statetxt = 'groups sorted by cell count, ascending';
            otherwise
                idx = (1:numel(rowtotal))';
                statetxt = '';
        end
    end

function in_callback_sortgroups(~, ~)
        % Descending, then ascending, then back to the order it started in.
        % Both tabs move together, each sorting its own groups by its own
        % counts, so the button means the same thing whichever tab is up.
        sortstate = mod(sortstate+1, 3);
        for kloop = 1:2
            in_drawtab(kloop);
        end
    end

function in_plot1(kp1, idxp1, stxtp1)
        axh = ax{kp1,1};
        y = Tall{kp1}(idxp1, :);
        cla(axh);
        b = in_stackedbar(axh, y);
        xticks(axh, 1:size(y, 1));
        labelsx1 = gui.i_escapeunderscore(labelsxall{kp1}(idxp1));
        xticklabels(axh, labelsx1);

        xlabel(axh, strrep(clabelall{kp1}, '_', '\_'))
        ylabel(axh, '# of cells')
        subtitle(axh, stxtp1);
        labelsy1 = gui.i_escapeunderscore(labelsyall{kp1});
        lgd = legend(axh, b, labelsy1, 'Location', 'bestoutside');
        title(lgd, strrep(llabelall{kp1}, '_', '\_'));
    end

function in_plot2(kp2, idxp2)
        axh = ax{kp2,2};
        Tk = Tall{kp2}(idxp2, :);
        cla(axh);
        labelsx2 = gui.i_escapeunderscore(labelsxall{kp2}(idxp2));
        clab = strrep(clabelall{kp2}, '_', '\_');
        if ismosaic
            hleg = in_drawmosaic(axh, Tk, labelsx2);
            xlabel(axh, [clab '  (bar width = group size)'])
        else
            hleg = in_drawpctbars(axh, Tk, labelsx2);
            xlabel(axh, clab)
        end
        ylabel(axh, '% of cells')
        ylim(axh, [0, 1]);
        if isempty(hleg), return; end
        labelsy2 = gui.i_escapeunderscore(labelsyall{kp2});
        lgd = legend(axh, hleg, labelsy2, 'Location', 'bestoutside');
        title(lgd, strrep(llabelall{kp2}, '_', '\_'));
    end

function hleg = in_drawpctbars(axh, Tk, labs)
        % Every group the same width, whatever its size.
        y = Tk ./ sum(Tk, 2);
        b = in_stackedbar(axh, y);
        xticks(axh, 1:size(y, 1));
        xticklabels(axh, labs);
        hleg = b;
    end

function b = in_stackedbar(axh, y)
        % Colour the series outright instead of through the colormap.
        % FaceColor 'flat' with a numeric CData resolves a colour from the
        % axes Colormap and CLim, which are per-axes: the counts panel and
        % the percentage panel are separate axes, and when their state drifts
        % apart the same category comes out a different colour in each.
        b = bar(axh, y, 'stacked');
        cols = in_seriescolors(size(y, 2));
        for kx = 1:size(y, 2)
            b(kx).FaceColor = cols(kx,:);
        end
    end

function cols = in_seriescolors(n)
        % One row per category, shared by every panel so a category keeps
        % its colour across the counts bars, the percentage bars and the
        % mosaic.
        cols = turbo(max(n, 2));
        cols = cols(1:max(n, 1), :);
    end

function hleg = in_drawmosaic(axh, Tk, labs)
        % Same stacked shares as IN_DRAWPCTBARS, but each column is as wide
        % as its group is big, so a group of three stops occupying as much
        % of the plot as a group of three thousand. BAR cannot vary width
        % per bar, so the columns are drawn as patches.
        rowtot = sum(Tk, 2);
        keep = rowtot > 0;
        hleg = gobjects(0);
        if ~any(keep), return; end

        m = size(Tk, 2);
        cols = in_seriescolors(m);
        gap = 0.006;
        usable = 1 - gap*(sum(keep) - 1);
        w = usable * rowtot / sum(rowtot);

        hleg = gobjects(m, 1);
        xc = nan(numel(rowtot), 1);
        x0 = 0;
        washeld = ishold(axh);
        hold(axh, 'on');
        for i = 1:numel(rowtot)
            if ~keep(i), continue; end
            p = Tk(i,:) / rowtot(i);
            yb = 0;
            for j = 1:m
                yt = yb + p(j);
                hp = patch(axh, [x0 x0+w(i) x0+w(i) x0], [yb yb yt yt], cols(j,:));
                if ~isgraphics(hleg(j)), hleg(j) = hp; end
                yb = yt;
            end
            xc(i) = x0 + w(i)/2;
            x0 = x0 + w(i) + gap;
        end
        if ~washeld, hold(axh, 'off'); end

        xlim(axh, [0, x0-gap]);
        xticks(axh, xc(keep));
        xticklabels(axh, labs(keep));
    end

function in_callback_togglemosaic(~, ~)
        % Applies to both tabs, like the sort, so the button means the same
        % thing whichever tab is up.
        ismosaic = ~ismosaic;
        for kmos = 1:2
            in_drawtab(kmos);
        end
    end

function in_callback_saveCrossTable(~, ~)
        % Save what is on screen: the tab being looked at, in the order it
        % is currently sorted, so the file matches the plot.
        ksave = 1;
        for ktsave = 1:2
            if isequal(tabgp.SelectedTab, tab{ktsave}), ksave = ktsave; end
        end
        idxsave = in_grouporder(ksave);
        t = array2table(Tall{ksave}(idxsave, :));
        try
            t.Properties.VariableNames = labelsyall{ksave};
            t.Properties.RowNames = labelsxall{ksave}(idxsave);
        catch
            % keep default headers if labels don't match table dimensions
        end
        gui.i_exporttable(t, true, 'Tcrosstabul', ...
            'CrosstabulTable',[],[], hx.FigHandle);
    end

function tf = i_iscontinuous(v)
        tf = isnumeric(v) && ~islogical(v) && ...
            numel(unique(v)) > min(20, numel(v) * 0.05);
    end
end
