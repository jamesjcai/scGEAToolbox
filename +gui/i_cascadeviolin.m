function i_cascadeviolin(sce, Xt, thisc, glist, ytxt, grouporder, colorit)

if nargin < 7, colorit = false; end
if nargin < 6, grouporder = []; end
if nargin < 5, ytxt = ''; end

F = cell(length(glist), 1);
for k = 1:length(glist)
    [~, idx] = ismember(glist(k), sce.g);
    y = full(Xt(idx, :));
    ttxt = sce.g(idx);

    f = figure('visible', 'off');
    ax1 = subplot(1, 2, 1);
    pkg.i_violinplot(y, thisc, colorit, grouporder);
    assignin('base', 'y', y);
    assignin('base', 'c', thisc);

    title(strrep(ttxt, '_', '\_'));
    ylabel(ytxt);
    subplot(1, 2, 2)
    sc_scattermarker(sce.X, sce.g, ...
        sce.s, glist(k), 2);


    tb = uitoolbar(f);
    pkg.i_addbutton2fig(tb, 'off', {@in_savedata, y, thisc}, ...
        'export.gif', 'Export data...');
    % pkg.i_addbutton2fig(tb,'off',{@gui.i_savemainfig,3}, ...
    %     "powerpoint.gif",'Save Figure to PowerPoint File...');
    pkg.i_addbutton2fig(tb, 'off', {@i_invertcolor, ax1, colorit, y, thisc, grouporder}, ...
        "xpowerpoint.gif", 'Switch to BW');
    P = get(f, 'Position');
    set(f, 'Position', [P(1) - 20 * k, P(2) - 20 * k, P(3), P(4)]);
    set(f, 'visible', 'on');
    f.Position(3) = f.Position(3) * 2.2;
    drawnow;
    F{k} = f;
end
gui.i_export2pptx(F, glist);


    function i_invertcolor(~, ~, ax1, colorit, y, thisc, grouporder)
        colorit = ~colorit;
        % delete(vh);
        axes(ax1)
        cla
        pkg.i_violinplot(y, thisc, colorit, grouporder);
end

end

function in_savedata(~, ~, a, b)
    T = table(a(:), b(:));
    T.Properties.VariableNames = {'ExprLevel', 'GroupID'};
    T = sortrows(T, 'ExprLevel', 'descend');
    T = sortrows(T, 'GroupID');
    gui.i_exporttable(T, true, 'Tviolindata','ViolinPlotTable');
end

% "Tcellattrib","CellAttribTable"
% "Tviolindata","ViolinPlotTable"