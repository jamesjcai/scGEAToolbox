function [hFig] = i_spiderplot(Y, thisc, labelx, sce, parentfig, levelorder, opts)
% LEVELORDER, optional, is the order the groups were listed in when the
% user picked them; the polygons and legend keep it. Groups it leaves out
% follow in sorted order.
% ValueKind="score" (default) or "expression" names what Y holds, in the
% save button, its dialog and the saved tables' default names.

arguments
    Y
    thisc
    labelx
    sce = []
    parentfig = []
    levelorder = []
    opts.ValueKind (1,1) string {mustBeMember(opts.ValueKind, ["score", "expression"])} = "score"
end

switch opts.ValueKind
    case "expression"
        saveTip = 'Save Expression Matrix...';
        saveWhat = 'expression';
        perCellOption = 'Per-cell Expression';
        perCellNames = {'Tgeneexpr', 'GeneExprTable'};
        meanNames = {'Tgeneexprmean', 'GeneExprMeanTable'};
    otherwise
        % "score", the only other value the arguments block admits.
        saveTip = 'Save Score Matrix...';
        saveWhat = 'score';
        perCellOption = 'Per-cell Scores';
        perCellNames = {'Tcellsignmt', 'CellSignatTable'};
        meanNames = {'Tspiderdata', 'SpiderOutTable'};
end
saveTitle = saveTip(1:end-3);

[c, cL] = findgroups(string(thisc));
if ~isempty(levelorder)
    [~, neworder] = ismember(string(levelorder), cL);
    neworder = neworder(neworder > 0);
    neworder = [neworder(:); setdiff((1:numel(cL)).', neworder(:))];
    [~, rank] = sort(neworder);
    c = rank(c);
    cL = cL(neworder);
end
% P = grpstats(Y, c, 'mean');
% Grouped means: rows = cell groups (spider series), columns = signatures
% (spider axes). Do not transpose -- spider_plot expects series-by-axes, and
% the axis labels (labelx) index the columns.
P = splitapply(@mean, Y, c);
n = size(P, 2);
% Keep the names as given for the exported table: the ones drawn below lose
% their "(Collection)" prefix and get TeX-escaped underscores.
rawlabels = cellstr(string(labelx));

%         axes_limits=[repmat(min([0, min(P(:))]),1,n);...
%             repmat(max(P(:)),1,n)];

axes_limits = [repmat(min(P(:)), 1, n); ...
repmat(max(P(:)), 1, n)];


if ~isempty(strfind(labelx{1}, ')'))
    titlex = extractBefore(labelx{1}, strfind(labelx{1}, ')')+1);
    labelx = extractAfter(labelx, strfind(labelx{1}, ')')+1);
else
    titlex = '';
end

hx=gui.myFigure(parentfig);
hFig=hx.FigHandle;
% Every draw below, the button callbacks included, targets this axes: the
% current axes after a dialog can be another window's.
ax = hx.AxHandle;

hx.addCustomButton('off', {@i_savedata}, 'floppy-disk-arrow-in.jpg', saveTip);
hx.addCustomButton('off', @i_showvalues, "heap_snapshot_large_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", 'Switch Values ON/OFF');
hx.addCustomButton('off', @i_reordersamples, "keyframe-plus-in.jpg", 'Switch Legend ON/OFF');
hx.addCustomButton('on', @i_editgrpnames, 'edit.jpg', 'Rename Group Names');

showaxes = true;
showlegend = true;

bkcolor = gui.i_getthemebkgcolor(hFig);

labelx = gui.i_escapeunderscore(labelx);
spider_plot_R2019b(P, 'AxesHandle', ax, 'AxesLabels', labelx, ...
'AxesPrecision', 2, 'AxesLimits', axes_limits, ...
'BackgroundColor', bkcolor, ...
'AxesFontColor', 1-bkcolor, ...
'AxesZeroColor', 1-bkcolor);
cL = gui.i_escapeunderscore(cL);
legend(ax, cL, 'Location', 'best');
if ~isempty(titlex), title(ax, titlex); end
hx.show(parentfig);


function i_showvalues(~, ~)
        showaxes = ~showaxes;
        cla(ax, 'reset');
        if showaxes
            spider_plot_R2019b(P, 'AxesHandle', ax, 'AxesLabels', labelx, ...
                'AxesDisplay', 'all', 'AxesPrecision', 2, ...
                'AxesLimits', axes_limits, ...
                'BackgroundColor', bkcolor, ...
                'AxesFontColor', 1-bkcolor, ...
                'AxesZeroColor', 1-bkcolor);
        else
            spider_plot_R2019b(P, 'AxesHandle', ax, 'AxesLabels', labelx, ...
                'AxesDisplay', 'none', 'AxesPrecision', 2, ...
                'AxesLimits', axes_limits,...
                'BackgroundColor', bkcolor, ...
                'AxesFontColor', 1-bkcolor, ...
                'AxesZeroColor', 1-bkcolor);
        end
        if showlegend, legend(ax, cL); end
        if ~isempty(titlex), title(ax, titlex); end
    end

function i_reordersamples(~, ~)
        showlegend = ~showlegend;

        if showlegend
            legend(ax, cL, 'Location', 'best');
        else
            legend(ax, 'off');
        end
    end


function i_editgrpnames(~, ~)

        if gui.i_isuifig(parentfig)
            [indxx, tfx] = gui.myListdlg(hFig, string(cL), 'Select group name', [], false);
        else
            [indxx, tfx] = listdlg('PromptString', ...
                {'Select group name'}, ...
                'SelectionMode', 'single', ...
                'ListString', string(cL), 'ListSize', [220, 300]);
        end

        if tfx == 1
            i = ismember(c, indxx);
            if gui.i_isuifig(parentfig)
                newctype = gui.myInputdlg({'New cell type'}, 'Rename', cL(c(i)), hFig);
            else
                newctype = inputdlg('New cell type', 'Rename', [1, 50], cL(c(i)));
            end
            if ~isempty(newctype)
                cL(c(i)) = newctype;
                legend(ax, cL, 'Location', 'best');
            end
        end
    end


function i_savedata(~, ~)
        answerx = gui.myQuestdlg(hFig, ['Save the ', saveWhat, ' of every cell, ' ...
            'or the group means the radar plot draws?'], saveTitle, ...
            {perCellOption, 'Group Means', 'Cancel'}, perCellOption);
        switch answerx
            case perCellOption
                if isempty(sce)
                    a = string(1:size(Y, 1));
                else
                    a = matlab.lang.makeUniqueStrings(string(sce.c_cell_id));
                end
                T = array2table(Y, 'VariableNames', rawlabels, 'RowNames', a);
                T.Cell_Group = strrep(string(cL(c)), '\_', '_');
                T.Properties.DimensionNames{1} = 'Cell_ID';
                gui.i_exporttable(T, false, ...
                    perCellNames{:}, [], [], hFig);
            case 'Group Means'
                % cL is TeX-escaped for the legend (and may have been renamed).
                grp = matlab.lang.makeUniqueStrings(strrep(string(cL), '\_', '_'));
                T = array2table(P, 'VariableNames', rawlabels, 'RowNames', grp);
                T.NumCells = accumarray(c(:), 1);
                T.Properties.DimensionNames{1} = 'Cell_Group';
                gui.i_exporttable(T, false, ...
                    meanNames{:}, [], [], hFig);
            otherwise
                return;
        end
    end

end
