function [hFig] = i_violinplot(y, thisc, ttxt, colorit, cLorder, posg, parentfig)
% see also: gui.sc_uitabgrpfig_violin

if nargin < 7, parentfig = []; end
if nargin < 6, posg = []; end
if nargin < 5 || isempty(cLorder), [~, cLorder] = grp2idx(thisc); end
if nargin < 4, colorit = true; end
if nargin < 3, ttxt = ''; end

if matlab.ui.internal.isUIFigure(parentfig), focus(parentfig); end
hFig = figure('visible', 'off', 'DockControls', 'off');

isdescend = false;
OldTitle = [];
% OldXTickLabel = [];
cLorder = strrep(cLorder, '_', '\_');
thisc = strrep(string(thisc), '_', '\_');
pkg.i_violinplot(y, thisc, colorit, cLorder);
title(strrep(ttxt, '_', '\_'));
%ylabel(selitems{indx1});


% tb = uitoolbar(hFig);
tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb, 'Separator', 'off');

pkg.i_addbutton2fig(tb, 'off', @i_savedata, ...
    'export.gif', 'Export data...');
pkg.i_addbutton2fig(tb, 'off', @i_testdata, ...
    'icon-fa-stack-exchange-10.gif', 'ANOVA/T-test...');
pkg.i_addbutton2fig(tb, 'off', @i_addsamplesize, ...
    "icon-mat-blur-linear-10.gif", 'Add Sample Size');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, ...
    "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'off', @i_invertcolor, ...
    "plotpicker-pie.gif", 'Switch BW/Color');
pkg.i_addbutton2fig(tb, 'off', @i_reordersamples, ...
    "plotpicker-errorbar.gif", 'Reorder Samples');

pkg.i_addbutton2fig(tb, 'off', @i_selectsamples, ...
    "plotpicker-errorbarx.gif", 'Select Samples');


pkg.i_addbutton2fig(tb, 'off', @i_sortbymean, ...
    "plotpicker-cra.gif", 'Sort Samples by Median');
pkg.i_addbutton2fig(tb, 'off', @i_renametitle, ...
    "icon-mat-touch-app-10.gif", 'Change Plot Title');
pkg.i_addbutton2fig(tb, 'on', @i_viewgenenames, ...
    'HDF_point.gif', 'Show Gene Names');
pkg.i_addbutton2fig(tb, 'on', @i_resizewin, ...
    'HDF_pointx.gif', 'Resize Plot Window');


i_addsamplesize([],[]);
i_testdata([],[]);
if nargout > 0
    return;
end

try
if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
    [px_new] = gui.i_getchildpos(parentfig, hFig);
    if ~isempty(px_new)
        movegui(hFig, px_new);
    else
        movegui(hFig, 'center');
    end
else
    movegui(hFig, 'center');
end
catch
    movegui(hFig, 'center');
end
set(hFig, 'visible', 'on');

%catch ME
%    errordlg(ME.message);
%end


    function i_resizewin(~,~)
        %oldw
        %oldh
        w = gui.i_inputnumk(450, 10, 2000, 'Window width');
        if isempty(w), return; end
        h = gui.i_inputnumk(420, 10, 2000, 'Window height');
        if isempty(h), return; end
        hFig.Position = [hFig.Position(1) hFig.Position(2) w h];
    end

    function i_renametitle(~, ~)
        helpdlg('Double-click on the title to make change.', '');
    end

    function i_invertcolor(~, ~)
        colorit = ~colorit;
        b = hFig.get("CurrentAxes");
        cla(b);
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end

    function i_addsamplesize(~, ~)
        % b = gca;
        b = hFig.get("CurrentAxes");
        b.FontName='Palatino';
        % assert(isequal(cLorder, b.XTickLabel));

        if isequal(cLorder, b.XTickLabel)
            a = zeros(length(cLorder), 1);            
            for k = 1:length(cLorder)
                a(k) = sum(thisc == cLorder(k));
                cb=pad([string(b.XTickLabel{k}); sprintf("(n=%d)",a(k))],'both');
                b.XTickLabel{k} = sprintf('%s\\newline%s', cb(:));                
            end
        else
            b.XTickLabel = cLorder;                
        end        
    end

    function i_sortbymean(~, ~)
        [cx, cLx] = grp2idx(thisc);
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

        %[~,cL,noanswer]=gui.i_reordergroups(thisc, cLx_sorted, f);
        %if noanswer, return; end
        b = hFig.get("CurrentAxes");
        cla(b);
        cLorder = cLx_sorted;
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_reordersamples(~, ~)
        [~, cLorder, noanswer] = gui.i_reordergroups(thisc, [], hFig);

        % cLorder
        if noanswer, return; end
        b = hFig.get("CurrentAxes");
        cla(b);
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_selectsamples(~, ~)
        [~,cL]=grp2idx(thisc);
        [newidx] = gui.i_selmultidlg(cL, cLorder);
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
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_testdata(~, ~)
        a = hFig.get("CurrentAxes");
        if isempty(OldTitle)
            OldTitle = a.Title.String;
            if size(y, 2) ~= length(thisc)
                y = y.';
            end
            tbl = pkg.e_grptest(y, thisc);
            %h1=gca;
            %titre=string(h1.Title.String);

            %     a=sprintf('%s\n%s=%.2e; %s=%.2e', ...
            %         strrep(string(ttxt),'_','\_'), ...
            %         strrep(tbl.Properties.VariableNames{1},'_','\_'), ...
            %         tbl.(tbl.Properties.VariableNames{1}), ...
            %         strrep(tbl.Properties.VariableNames{2},'_','\_'), ...
            %         tbl.(tbl.Properties.VariableNames{2}));
            %     title(a);
            if ~isempty(tbl) && istable(tbl)
                b = sprintf('%s=%.2e; %s=%.2e', ...
                    strrep(tbl.Properties.VariableNames{1}, '_', '\_'), ...
                    tbl.(tbl.Properties.VariableNames{1}), ...
                    strrep(tbl.Properties.VariableNames{2}, '_', '\_'), ...
                    tbl.(tbl.Properties.VariableNames{2}));
            else
                % b='p\_ttest=N.A.; p\_wilcoxon=N.A.';
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
            helpdlg('The gene set is empty. This score may not be associated with any gene set.');
        else
            %idx=matches(sce.g,posg,'IgnoreCase',true);
            %gg=sce.g(idx);
            inputdlg(ttxt, ...
                '', [10, 50], ...
                {char(posg)});
        end
    end

    function i_savedata(~, ~)    
        T = table(y(:), thisc(:));
        T.Properties.VariableNames = {'ScoreLevel', 'GroupID'};
        %T=sortrows(T,'ScoreLevel','descend');
        %T=sortrows(T,'GroupID');
        gui.i_exporttable(T, true, 'Tviolindata','ViolinPlotTable');
    end

end
