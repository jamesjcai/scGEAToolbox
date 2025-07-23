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

hx=gui.myFigure;
hFig=hx.FigHandle;

hx.addCustomButton('off', {@i_savedata}, 'floppy-disk-arrow-in.jpg', 'Export data...');
hx.addCustomButton('off', @i_showvalues, "heap_snapshot_large_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", 'Switch Values ON/OFF');
hx.addCustomButton('off', @i_reordersamples, "keyframe-plus-in.jpg", 'Switch Legend ON/OFF');
hx.addCustomButton('on', @i_editgrpnames, 'edit.jpg', 'Rename Group Names');

showaxes = true;
showlegend = true;

bkcolor = gui.i_getthemebkgcolor(hFig);

labelx = strrep(labelx, '_', '\_');
spider_plot_R2019b(P, 'AxesLabels', labelx, ...
    'AxesPrecision', 2, 'AxesLimits', axes_limits, ...
    'BackgroundColor', bkcolor, ...
    'AxesFontColor', 1-bkcolor, ...
    'AxesZeroColor', 1-bkcolor);
cL = strrep(cL, '_', '\_');
legend(cL, 'Location', 'best');
if ~isempty(titlex), title(titlex); end

hx.show(parentfig);


%catch ME
%    gui.myErrordlg(parentfig, ME.message, ME.identifier);
%end


    function i_showvalues(~, ~)
        showaxes = ~showaxes;
        cla(hx.AxHandle);
        if showaxes
            spider_plot_R2019b(P, 'AxesLabels', labelx, ...
                'AxesDisplay', 'all', 'AxesPrecision', 2, ...
                'AxesLimits', axes_limits, ...
                'BackgroundColor', bkcolor, ...
                'AxesFontColor', 1-bkcolor, ...
                'AxesZeroColor', 1-bkcolor);
        else
            spider_plot_R2019b(P, 'AxesLabels', labelx, ...
                'AxesDisplay', 'none', 'AxesPrecision', 2, ...
                'AxesLimits', axes_limits,...
                'BackgroundColor', bkcolor, ...
                'AxesFontColor', 1-bkcolor, ...
                'AxesZeroColor', 1-bkcolor);
        end
        if showlegend, legend(cL); end
        if ~isempty(titlex), title(titlex); end
    end

    % function i_addsamplesize(~, ~)
    %     % b=gca;
    %     % if isempty(OldXTickLabel)
    %     %     a=zeros(length(cLorder),1);
    %     %
    %     %     OldXTickLabel=b.XTickLabel;
    %     %     for k=1:length(cLorder)
    %     %         a(k)=sum(thisc==cLorder(k));
    %     %         b.XTickLabel{k}=sprintf('%s\\newline(n=%d)', ...
    %     %             b.XTickLabel{k},a(k));
    %     %     end
    %     % else
    %     %     b.XTickLabel=OldXTickLabel;
    %     %     OldXTickLabel=[];
    %     % end
    % end

%    function i_sortbymean(~, ~)
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
%end


    function i_reordersamples(~, ~)
        showlegend = ~showlegend;

        if showlegend
            legend(cL, 'Location', 'best');
        else
            legend off
        end
    end


    function i_editgrpnames(~, ~)

        if gui.i_isuifig(parentfig)
            [indxx, tfx] = gui.myListdlg(parentfig, string(cL), 'Select group name');
        else
            [indxx, tfx] = listdlg('PromptString', ...
                {'Select group name'}, ...
                'SelectionMode', 'single', ...
                'ListString', string(cL), 'ListSize', [220, 300]);
        end

        if tfx == 1
            i = ismember(c, indxx);
            if gui.i_isuifig(parentfig)
                newctype = gui.myInputdlg({'New cell type'}, 'Rename', cL(c(i)), parentfig);
            else
                newctype = inputdlg('New cell type', 'Rename', [1, 50], cL(c(i)));
            end
            if ~isempty(newctype)
                cL(c(i)) = newctype;
                legend(cL, 'Location', 'best');
            end
        end
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
