function [hFig] = i_spiderplot(Y, thisc, labelx, sce, parentfig)

if nargin < 5, parentfig = []; end
if nargin < 4, sce = []; end

[c, cL] = grp2idx(thisc);
P = grpstats(Y, c, 'mean');
n = size(P, 2);

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

hFig = figure('visible', 'off');
% tb = uitoolbar(hFig);

tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb, 'Separator', 'off');

pkg.i_addbutton2fig(tb, 'off', {@i_savedata}, ...
    'export.gif', 'Export data...');
% pkg.i_addbutton2fig(tb,'off',{@i_testdata,P,thisc}, ...
%     'icon-fa-stack-exchange-10.gif','ANOVA/T-test...');
% pkg.i_addbutton2fig(tb,'off',@i_addsamplesize, ...
%     "icon-mat-blur-linear-10.gif",'Add Sample Size');
pkg.i_addbutton2fig(tb, 'on', {@gui.i_savemainfig, 3}, ...
    "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'off', @i_invertcolor, ...
    "plotpicker-pie.gif", 'Switch Values ON/OFF');
pkg.i_addbutton2fig(tb, 'off', @i_reordersamples, ...
    "plotpicker-errorbar.gif", 'Switch Legend ON/OFF');
% pkg.i_addbutton2fig(tb,'off',@i_sortbymean, ...
%     "plotpicker-cra.gif",'Sort Samples by Median');
pkg.i_addbutton2fig(tb, 'off', @i_renametitle, ...
    "icon-mat-touch-app-10.gif", 'Change Plot Title');
pkg.i_addbutton2fig(tb, 'on', @i_viewgenenames, ...
    'HDF_point.gif', 'Rename Group Names');

showaxes = true;
showlegend = true;

labelx = strrep(labelx, '_', '\_');
spider_plot_R2019b(P, 'AxesLabels', labelx, ...
    'AxesPrecision', 2, 'AxesLimits', axes_limits);
cL = strrep(cL, '_', '\_');
legend(cL, 'Location', 'best');
if ~isempty(titlex), title(titlex); end

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
set(hFig, 'visible', 'on');

%catch ME
%    errordlg(ME.message);
%end

    function i_renametitle(~, ~)
        helpdlg('Double-click on the title to make change.', '');
    end

    function i_invertcolor(~, ~)
        showaxes = ~showaxes;
        cla;
        if showaxes
            spider_plot_R2019b(P, 'AxesLabels', labelx, ...
                'AxesDisplay', 'all', 'AxesPrecision', 2, ...
                'AxesLimits', axes_limits);
        else
            spider_plot_R2019b(P, 'AxesLabels', labelx, ...
                'AxesDisplay', 'none', 'AxesPrecision', 2, ...
                'AxesLimits', axes_limits);
        end
        if showlegend, legend(cL); end
        if ~isempty(titlex), title(titlex); end
    end

    function i_addsamplesize(~, ~)
        % b=gca;
        % if isempty(OldXTickLabel)
        %     a=zeros(length(cLorder),1);
        %
        %     OldXTickLabel=b.XTickLabel;
        %     for k=1:length(cLorder)
        %         a(k)=sum(thisc==cLorder(k));
        %         b.XTickLabel{k}=sprintf('%s\\newline(n=%d)', ...
        %             b.XTickLabel{k},a(k));
        %     end
        % else
        %     b.XTickLabel=OldXTickLabel;
        %     OldXTickLabel=[];
        % end
    end

    function i_sortbymean(~, ~)
        % [cx,cLx]=grp2idx(thisc);
        % a=zeros(max(cx),1);
        % for k=1:max(cx)
        %     a(k)=median(y(cx==k));
        % end
        % if isdescend
        %     [~,idx]=sort(a,'ascend');
        %     isdescend=false;
        % else
        %     [~,idx]=sort(a,'descend');
        %     isdescend=true;
        % end
        % cLx_sorted=cLx(idx);
        %
        % %[~,cL,noanswer]=gui.i_reordergroups(thisc,cLx_sorted);
        % %if noanswer, return; end
        % cla
        % cLorder=cLx_sorted;
        % pkg.i_violinplot(y,thisc,colorit,cLorder);
end


    function i_reordersamples(~, ~)
        showlegend = ~showlegend;

        if showlegend
            legend(cL, 'Location', 'best');
        else
            legend off
        end
    end


    function i_viewgenenames(~, ~)
        [indxx, tfx] = listdlg('PromptString', ...
            {'Select group name'}, ...
            'SelectionMode', 'single', ...
            'ListString', string(cL));
        if tfx == 1
            i = ismember(c, indxx);
            newctype = inputdlg('New cell type', 'Rename', [1, 50], cL(c(i)));
            if ~isempty(newctype)
                cL(c(i)) = newctype;
                legend(cL, 'Location', 'best');
            end
        end
        % if isempty(posg)
        %     helpdlg('The gene set is empty. This score may not be associated with any gene set.');
        % else
        %     %idx=matches(sce.g,posg,'IgnoreCase',true);
        %     %gg=sce.g(idx);
        %     inputdlg(ttxt, ...
        %         '',[10 50], ...
        %         {char(posg)});
        % end
    end


    function i_savedata(~, ~)
        if isempty(sce)
            a = string(1:size(Y, 1));
        else
            a = matlab.lang.makeUniqueStrings(sce.c_cell_id);
        end
        T = array2table(Y, 'VariableNames', ...
            labelx, 'RowNames', a);
        name = 'Cell_Group';
        %T.(name) = thisc;
        T.(name) = cL(c);
        T.Properties.DimensionNames{1} = 'Cell_ID';
        needwait = false;                                    
        gui.i_exporttable(T, needwait, ...
            'Tspiderdata', 'SpiderOutTable');
    end

end
