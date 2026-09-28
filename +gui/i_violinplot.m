function [hFig] = i_violinplot(y, thisc, ttxt, colorit, cLorder, ...
    posg, parentfig)


if nargin < 7, parentfig = []; end
if nargin < 6, posg = []; end
if nargin < 5 || isempty(cLorder), [~, cLorder] = findgroups(string(thisc)); end
if nargin < 4, colorit = true; end
if nargin < 3, ttxt = ''; end

if gui.i_isuifig(parentfig), focus(parentfig); end
hx=gui.myFigure(parentfig);
hFig=hx.FigHandle;


isdescend = false;
OldTitle = [];
cLorder = gui.i_escapeunderscore(cLorder);
thisc = strrep(string(thisc), '_', '\_');
pkg.i_bindviolinplot(y, thisc, colorit, cLorder, hx.AxHandle);
title(hx.AxHandle, gui.i_escapeunderscore(ttxt));


hx.addCustomButton('off', @i_savedata, 'floppy-disk-arrow-in.jpg', 'Export data...');
hx.addCustomButton('off', @i_testdata, 'mw-pickaxe-mining.jpg', 'ANOVA/T-test...');
hx.addCustomButton('off', @i_addsamplesize, "unjoin3d.jpg", 'Add Sample Size');
hx.addCustomButton('off', @i_invertcolor, "align-top-box-solid.jpg", 'Switch BW/Color');
hx.addCustomButton('off', @i_reordersamples, "rebase_edit_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", 'Reorder Samples');
hx.addCustomButton('off', @i_selectsamples, "edit.jpg", 'Select Samples');
hx.addCustomButton('off', @i_sortbymean, "reorder.jpg", 'Sort Samples by Median');
hx.addCustomButton('on', @i_viewgenenames, 'HDF_point.gif', 'Show Gene Names');


if nargout > 0, return; end
hx.show(parentfig);


function i_invertcolor(~, ~)
        colorit = ~colorit;
        b = hFig.get("CurrentAxes");
        cla(b);
        pkg.i_bindviolinplot(y, thisc, colorit, cLorder, b);
    end

function i_addsamplesize(~, ~)
        b = hFig.get("CurrentAxes");
        b.FontName='Palatino';

        if isequal(cLorder, b.XTickLabel)
            a = zeros(length(cLorder), 1);
            for k = 1:length(cLorder)
                a(k) = sum(thisc == cLorder(k));
                cb=pad([string(b.XTickLabel{k}); sprintf("(n=% d)",a(k))],'both');
                b.XTickLabel{k} = sprintf('%s\\newline%s', cb(:));
            end
        else
            b.XTickLabel = cLorder;
        end
    end

function i_sortbymean(~, ~)
        [cx, cLx] = findgroups(string(thisc));
        a = zeros(max(cx), 1);
        for k = 1:max(cx)
            a(k) = median(y(cx == k));
        end
        if isdescend
            [~, idx] = sort(a, 'ascend');
            isdescend = false;
        else
            [~, idx] = sort(a, 'descend');
            isdescend = true;
        end
        cLx_sorted = cLx(idx);

        b = hFig.get("CurrentAxes");
        cla(b);
        cLorder = cLx_sorted;
        pkg.i_bindviolinplot(y, thisc, colorit, cLorder, b);
    end


function i_reordersamples(~, ~)
        [~, cLorder, noanswer] = gui.i_reordergroups(thisc, [], hFig);

        % cLorder
        if noanswer, return; end
        b = hFig.get("CurrentAxes");
        cla(b);
        pkg.i_bindviolinplot(y, thisc, colorit, cLorder, b);
    end


function i_selectsamples(~, ~)
        [~,cL] = findgroups(string(thisc));
        % Parent the dialog on this plot window so it, not the main app, is
        % raised again when the dialog closes.
        [newidx] = gui.i_selmultidialog(cL, cLorder, hFig);
        if isempty(newidx), return; end
        picked=ismember(thisc,cL(newidx));
%        [~, cLorder, noanswer] = gui.i_reordergroups(thisc, [], f);
%        % cLorder
%        if noanswer, return; end

        cLorder=cLorder(ismember(cLorder,cL(newidx)));
        b = hFig.get("CurrentAxes");
        cla(b);
        y=y(picked);
        thisc=thisc(picked);
        pkg.i_bindviolinplot(y, thisc, colorit, cLorder, b);
    end


function i_testdata(~, ~)
        a = hFig.get("CurrentAxes");
        if isempty(OldTitle)
            OldTitle = a.Title.String;
            if size(y, 2) ~= length(thisc)
                y = y.';
            end
            tbl = pkg.e_grptest(y, thisc);

            if ~isempty(tbl) && istable(tbl)
                b = sprintf('%s=%.2e; %s=%.2e', ...
                    strrep(tbl.Properties.VariableNames{1}, '_', '\_'), ...
                    tbl.(tbl.Properties.VariableNames{1}), ...
                    strrep(tbl.Properties.VariableNames{2}, '_', '\_'), ...
                    tbl.(tbl.Properties.VariableNames{2}));
            else
                b='p_{ttest}=N.A.; p_{wilcoxon}=N.A.';
            end

            if iscell(OldTitle)
                newtitle = OldTitle;
            else
                newtitle = {OldTitle};
            end
            newtitle{2} = b;
            a.Title.String = newtitle;
        else
            a.Title.String = OldTitle;
            OldTitle = [];

        end
    end


function i_viewgenenames(~, ~)
        if isempty(posg)
            gui.myHelpdlg(hFig, ['The gene set is empty. This score ' ...
                'may not be associated with any gene set.']);
        else
            if gui.i_isuifig(parentfig)
                gui.myInputdlg({ttxt}, '', {char(posg)}, hFig);
            else
                inputdlg(ttxt, '', [15, 80], {char(posg)});
            end
        end
    end

function i_savedata(~, ~)
        T = table(y(:), thisc(:));
        T.Properties.VariableNames = {'ScoreLevel', 'GroupID'};
        gui.i_exporttable(T, true, 'Tviolindata', 'ViolinPlotTable', [], [], hFig);
    end

end
