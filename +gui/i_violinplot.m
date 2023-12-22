function [f] = i_violinplot(y, thisc, ttxt, colorit, cLorder, posg)

if nargin < 6, posg = []; end % used in callback_CompareGeneBtwCls
if nargin < 5 || isempty(cLorder)
    [~, cLorder] = grp2idx(thisc);
end
if nargin < 4
    colorit = true;
end
f = figure('visible', 'off');

isdescend = false;
OldTitle = [];
% OldXTickLabel = [];
cLorder = strrep(cLorder, '_', '\_');
thisc = strrep(string(thisc), '_', '\_');
pkg.i_violinplot(y, thisc, colorit, cLorder);
title(strrep(ttxt, '_', '\_'));
%ylabel(selitems{indx1});


tb = uitoolbar(f);
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



movegui(f, 'center');

i_addsamplesize([],[]);
i_testdata([],[]);
if nargout > 0
    return; 
end

set(f, 'visible', 'on');

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
        f.Position = [f.Position(1) f.Position(2) w h];
    end

    function i_renametitle(~, ~)
        helpdlg('Double-click on the title to make change.', '');
    end

    function i_invertcolor(~, ~)
        colorit = ~colorit;
        b = f.get("CurrentAxes");
        cla(b);
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end

    function i_addsamplesize(~, ~)
        % b = gca;
        b = f.get("CurrentAxes");
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

        %[~,cL,noanswer]=gui.i_reordergroups(thisc,cLx_sorted);
        %if noanswer, return; end
        b = f.get("CurrentAxes");
        cla(b);
        cLorder = cLx_sorted;
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_reordersamples(~, ~)
        [~, cLorder, noanswer] = gui.i_reordergroups(thisc);

        % cLorder
        if noanswer, return; end
        b = f.get("CurrentAxes");
        cla(b);
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_selectsamples(~, ~)
        [~,cL]=grp2idx(thisc);
        [newidx] = gui.i_selmultidlg(cL, cLorder);
        if isempty(newidx), return; end
        picked=ismember(thisc,cL(newidx));
%        [~, cLorder, noanswer] = gui.i_reordergroups(thisc);
%        % cLorder
%        if noanswer, return; end
        
        cLorder=cLorder(ismember(cLorder,cL(newidx)));
        b = f.get("CurrentAxes");
        cla(b);
        y=y(picked);
        thisc=thisc(picked);
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_testdata(~, ~)
        a = f.get("CurrentAxes");
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

            b = sprintf('%s=%.2e; %s=%.2e', ...
                strrep(tbl.Properties.VariableNames{1}, '_', '\_'), ...
                tbl.(tbl.Properties.VariableNames{1}), ...
                strrep(tbl.Properties.VariableNames{2}, '_', '\_'), ...
                tbl.(tbl.Properties.VariableNames{2}));
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