function sc_uitabgrpfig_vioplot(y, tabnamelist, thisc, parentfig)

if ~iscell(y), y = {y}; end
if nargin<4, parentfig = []; end
tabnamelist = string(tabnamelist);

[~, cLorder] = findgroups(string(thisc));
cLorder = strrep(cLorder, '_', '\_');

fw = gui.myWaitbar(parentfig);
isdescend = false;

thisc = strrep(string(thisc), '_', '\_');
colorit = true;

import mlreportgen.ppt.*;
% pw1 = fileparts(mfilename('fullpath'));
% pth = fullfile(pw1, '..', 'assets', 'Misc', 'myTemplate.pptx');

hx=gui.myFigure(parentfig);
hFig=hx.FigHandle;
hFig.Position(3) = hFig.Position(3) * 1.8;

%delete(findall(hFig, 'Tag', 'FigureToolBar'));
% if ~isempty(cx)
%     px = hFig.Position;
%     px_new = [cx(1)-px(3)/2 cx(2)-px(4)/2];
% else
%     px_new = [];
% end

n = length(tabnamelist);

tabgp = uitabgroup();

idx = 1;
focalg = tabnamelist(idx);
tab=cell(n,1);
ax0=cell(n,1);

OldTitle = cell(n,1);
for k=1:n
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',tabnamelist(k)));
    ax0{k} = axes('parent',tab{k});
    pkg.i_violinplot(y{k}, thisc, true, cLorder);
    title(ax0{k}, strrep(tabnamelist(k), '_', '\_'));
    % subtitle(ax0{k}, gui.i_getsubtitle(c));
    % gui.i_setautumncolor(c, a, true, any(c==0));
end
  

tabgp.SelectionChangedFcn = @displaySelection;


hx.addCustomButton('off',  @in_callback_genecards, 'www.jpg', 'GeneCards...');
hx.addCustomButton('off', @in_callback_showbarplot, "bar_chart_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", 'Switch to Bar Plot');
hx.addCustomButton('on', @in_callback_savedata, 'floppy-disk.jpg', 'Export summary data...');
hx.addCustomButton('off', @in_callback_savedata_alltab, 'floppy-disk-arrow-in.jpg', 'Export individual cell data... (new format)');
hx.addCustomButton('on', @in_callback_testdata, 'mw-pickaxe-mining.jpg', 'ANOVA/T-test...');
hx.addCustomButton('off', @in_callback_addsamplesize, "unjoin3d.jpg", 'Add Sample Size');
hx.addCustomButton('off', @in_callback_invertcolor, "align-top-box-solid.jpg", 'Switch BW/Color');
hx.addCustomButton('off', @in_callback_reordersamples, "rebase_edit_16dp_000000_FILL0_wght400_GRAD0_opsz20.jpg", 'Reorder Samples');
hx.addCustomButton('off', @in_callback_selectsamples, "edit.jpg", 'Select Samples');
hx.addCustomButton('off', @in_callback_sortbymean, "reorder.jpg", 'Sort Samples by Median');
hx.addCustomButton('off', @in_mergetabs, 'Brightness-3--Streamline-Core.jpg', 'Show on the same figure...');

hx.show(parentfig);
gui.myWaitbar(parentfig, fw);

ccx = true;


    function in_mergetabs(~, ~)
        figure;
        for kx = 1:n
            hAx2 = nexttile;
            hAx1 = ax0{kx};
            gui.i_cloneaxes(hAx1, hAx2);
        end
    end

    function in_callback_showbarplot(~,~)
        [~, idx]=ismember(focalg, tabnamelist); 
        [cx, cLx] = findgroups(string(thisc));
        a = zeros(max(cx), 1);
        for ks = 1:max(cx)
            a(ks) = median(y{idx}(cx == ks));
        end
        delete(ax0{idx});
        ax0{idx} = axes('parent',tab{idx});
        
        % assignin("base","y",y{idx});
        % assignin("base","thisc",thisc);

        %mv = grpstats(y{idx},thisc,@mean);
        mv = splitapply(@mean, y{idx}', cx);
        if ccx
            bar(mv,'w');            
        else
            bar(mv);
        end
        ccx = ~ccx;
        hold on
        %sv = splitapply(@std, y{idx}', cx)./sqrt(splitapply(@numel, y{idx}', cx));
        sv = grpstats(y{idx}, thisc, @std)./sqrt(grpstats(y{idx}, thisc, @numel));
        errorbar(1:length(mv), mv, zeros(size(sv)), sv, 'color', 'k' ,'linestyle','none');
        
        disp('Error bar shows the standard error of the mean (SEM), i.e., the standard deviation and dividing it by the square root of the sample size')
        set(ax0{idx},'xticklabel',cLx);
        
        title(ax0{idx}, strrep(tabnamelist(idx), '_', '\_'));
        if length(tab)>1
            answer = gui.myQuestdlg(hFig, 'Apply to other tabs?','');
            if ~strcmp(answer,'Yes'), return; end
            i_updatebarplot(idx);
        end
    end

    function i_updatebarplot(idx)
        if nargin<1, idx=[]; end
        [cx, cLx] = findgroups(string(thisc));
        for ks = 1:n
            if ks~=idx
                delete(ax0{ks});
                ax0{ks} = axes('parent',tab{ks});
                
                mv = grpstats(y{ks},thisc,@mean);
                % mv = splitapply(@mean, y{ks}, thisc);
                if ~ccx
                    colc = 'w';            
                else
                    colc = '';
                end
                bar(ax0{ks}, mv, colc);
                hold on
                sv = grpstats(y{ks}, thisc, @std)./sqrt(grpstats(y{ks}, thisc, @numel));
                % sv = splitapply(@std, y{ks}, thisc)./sqrt(splitapply(@numel, y{ks}, thisc));
                errorbar(ax0{ks}, 1:length(mv), mv, zeros(size(sv)), sv, 'color', 'k' ,'linestyle','none');
                set(ax0{ks},'xticklabel',cLx);
                title(ax0{ks}, strrep(tabnamelist(ks), '_', '\_'));  
            end
        end    
    end

    function displaySelection(~, event)
        t = event.NewValue;
        txt = t.Title;
        % disp("Viewing gene " + txt);
        [~,idx]=ismember(txt,tabnamelist);
        focalg = tabnamelist(idx);
    end

    function in_callback_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end

    function in_callback_updatealltab(idx)
        if nargin<1, idx = []; end
        for ks=1:n
            if ks~=idx
                delete(ax0{ks});
                ax0{ks} = axes('parent',tab{ks});
                pkg.i_violinplot(y{ks}, thisc, colorit, cLorder);
                title(ax0{ks}, strrep(tabnamelist(ks), '_', '\_'));           
            end
        end        
    end

    function in_callback_invertcolor(~, ~)
        colorit = ~colorit;
        [~,idx]=ismember(focalg, tabnamelist);
        delete(ax0{idx});
        ax0{idx} = axes('parent',tab{idx});
        pkg.i_violinplot(y{idx}, thisc, colorit, cLorder);
        title(ax0{idx}, strrep(tabnamelist(idx), '_', '\_'));
        tabgp.SelectedTab=tab{idx};
        drawnow;
        if length(tab)>1
            answer = gui.myQuestdlg(hFig, 'Apply to other tabs?','');
            if ~strcmp(answer,'Yes'), return; end
            in_callback_updatealltab(idx);
        end
    end

    function in_callback_updatesamplesizelabel(idx)
        if nargin<1, idx=[]; end
        for ks = 1:n
            if ks~=idx
                b = ax0{ks};        
                b.FontName='Palatino';
                if isequal(cLorder, b.XTickLabel)
                    a = zeros(length(cLorder), 1);            
                    for kx = 1:length(cLorder)
                        a(kx) = sum(thisc == cLorder(kx));
                        cb=pad([string(b.XTickLabel{kx}); sprintf("(n=%d)",a(kx))],'both');
                        b.XTickLabel{kx} = sprintf('%s\\newline%s', cb(:));                
                    end
                else
                    b.XTickLabel = cLorder;
                end
            end
        end
    end

    function in_callback_addsamplesize(~, ~)
        [~,idx]=ismember(focalg, tabnamelist);
        b = ax0{idx};
        b.FontName='Palatino';
        if isequal(cLorder, b.XTickLabel)
            a = zeros(length(cLorder), 1);            
            for kx = 1:length(cLorder)
                a(kx) = sum(thisc == cLorder(kx));
                cb=pad([string(b.XTickLabel{kx}); sprintf("(n=%d)",a(kx))],'both');
                b.XTickLabel{kx} = sprintf('%s\\newline%s', cb(:));            
            end
        else
            b.XTickLabel = cLorder;                
        end
        if length(tab)>1
            answer = gui.myQuestdlg(hFig, 'Apply to other tabs?','');
            if ~strcmp(answer,'Yes'), return; end
            in_callback_updatesamplesizelabel(idx);
        end
    end

    function in_callback_sortbymean(~, ~)
        [~,idx]=ismember(focalg, tabnamelist);       
        [cx, cLx] = findgroups(string(thisc));

        a = zeros(max(cx), 1);
        for ks = 1:max(cx)
            a(ks) = median(y{idx}(cx == ks));
        end
        if isdescend
            [~, idxx] = sort(a, 'ascend');
            isdescend = false;
        else
            [~, idxx] = sort(a, 'descend');
            isdescend = true;
        end
        cLx_sorted = cLx(idxx);
        
        if isequal(cLx, cLx_sorted)
           gui.myHelpdlg(hFig, 'Groups has already been sorted.');
        else
            delete(ax0{idx});
            ax0{idx} = axes('parent',tab{idx});       
            cLorder = cLx_sorted;
            pkg.i_violinplot(y{idx}, thisc, colorit, cLorder);
            title(ax0{idx}, strrep(tabnamelist(idx), '_', '\_'));
        end
    end

    function in_callback_reordersamples(~, ~)
        [~, cLorderx, noanswer] = gui.i_reordergroups(thisc);
        if noanswer, return; end
        [~,idx] = ismember(focalg, tabnamelist);
        delete(ax0{idx});
        ax0{idx} = axes('parent',tab{idx});
        pkg.i_violinplot(y{idx}, thisc, colorit, cLorderx);
        title(ax0{idx}, strrep(tabnamelist(idx), '_', '\_'));  

        if length(tab)>1
            answer = gui.myQuestdlg(hFig, 'Apply to other tabs?','');
            if ~strcmp(answer,'Yes'), return; end
            cLorder = cLorderx;
            in_callback_updatealltab(idx);
        end
    end

    function in_callback_selectsamples(~, ~)
        [~, cLorder] = findgroups(string(thisc));
        [newidx] = gui.i_selmultidialog(cLorder, cLorder, hFig);
        if isempty(newidx), return; end
        picked=ismember(thisc, cLorder(newidx));
       
        cLorderx = cLorder(ismember(cLorder,cLorder(newidx)));
        [~,idx]=ismember(focalg, tabnamelist);
        delete(ax0{idx});
        ax0{idx} = axes('parent',tab{idx});
        y_picked = y{idx}(picked);
        thisc_picked = thisc(picked);
        pkg.i_violinplot(y_picked, thisc_picked, colorit, cLorderx);
        title(ax0{idx}, strrep(tabnamelist(idx), '_', '\_'));

        if length(tab)>1        
            answer = gui.myQuestdlg(hFig, 'Apply to other tabs?','');
            if ~strcmp(answer,'Yes'), return; end
    
            for ks=1:n
                y{ks} = y{ks}(picked);
            end
            thisc = thisc_picked;
            cLorder = cLorderx;
            in_callback_updatealltab(idx);
        end
    end

    function in_callback_testdata(~, ~)
        for tabidx=1:n
            tabgp.SelectedTab=tab{tabidx};
            a = ax0{tabidx};
            thisy = y{tabidx};
            %a = hFig.get("CurrentAxes");
            if isempty(OldTitle{tabidx})
                OldTitle{tabidx} = a.Title.String;
                if size(thisy, 2) ~= length(thisc)
                    thisy = thisy.';
                end
                tbl = pkg.e_grptest(thisy, thisc);
                if ~isempty(tbl) && istable(tbl)
                    b = sprintf('%s = %.2e; %s = %.2e', ...
                        tbl.Properties.VariableNames{1}, ...
                        tbl.(tbl.Properties.VariableNames{1}), ...
                        tbl.Properties.VariableNames{2}, ...
                        tbl.(tbl.Properties.VariableNames{2}));
                else
                    if length(unique(thisc)) == 2
                        b='p_{ttest} = N.A.; p_{wilcoxon} = N.A.';
                    else
                        b='p_{anova} = N.A.; p_{kruskalwallis} = N.A.';
                    end
                end    
                if iscell(OldTitle{tabidx})
                    newtitle = OldTitle{tabidx};
                else
                    newtitle = OldTitle(tabidx);
                end
                newtitle{2} = b;
                a.Title.String = newtitle;
            else
                a.Title.String = OldTitle{tabidx};
                OldTitle{tabidx} = [];
            end
        end
        [~,tabidx]=ismember(focalg, tabnamelist);
        tabgp.SelectedTab=tab{tabidx};
    end
        
    % function i_testdataone(~,~)
    %     [~,idx]=ismember(focalg, tabnamelist);
    %     %delete(ax0{idx});
    %     %ax0{idx} = axes('parent',tab{idx});
    %     tabgp.SelectedTab=tab{idx};
    %     a = ax0{idx};
    %     thisy = y{idx};
    %     %a = hFig.get("CurrentAxes");
    %     if isempty(OldTitle{idx})
    %         OldTitle{idx} = a.Title.String;
    %         if size(thisy, 2) ~= length(thisc)
    %             thisy = thisy.';
    %         end
    %         tbl = pkg.e_grptest(thisy, thisc);
    %         %h1=gca;
    %         %titre=string(h1.Title.String);
    % 
    %         %     a=sprintf('%s\n%s=%.2e; %s=%.2e', ...
    %         %         strrep(string(ttxt),'_','\_'), ...
    %         %         strrep(tbl.Properties.VariableNames{1},'_','\_'), ...
    %         %         tbl.(tbl.Properties.VariableNames{1}), ...
    %         %         strrep(tbl.Properties.VariableNames{2},'_','\_'), ...
    %         %         tbl.(tbl.Properties.VariableNames{2}));
    %         %     title(a);
    %         if ~isempty(tbl) && istable(tbl)
    %             b = sprintf('%s=%.2e; %s=%.2e', ...
    %                 tbl.Properties.VariableNames{1}, ...
    %                 tbl.(tbl.Properties.VariableNames{1}), ...
    %                 tbl.Properties.VariableNames{2}, ...
    %                 tbl.(tbl.Properties.VariableNames{2}));
    % 
    %             % b = sprintf('%s=%.2e; %s=%.2e', ...
    %             %     strrep(tbl.Properties.VariableNames{1}, '_', '_'), ...
    %             %     tbl.(tbl.Properties.VariableNames{1}), ...
    %             %     strrep(tbl.Properties.VariableNames{2}, '_', '\_'), ...
    %             %     tbl.(tbl.Properties.VariableNames{2}));
    %         else
    %             if length(unique(thisc)) == 2
    %                 b='p_{ttest}=N.A.; p_{wilcoxon}=N.A.';
    %             else
    %                 b='p_{anova}=N.A.; p_{kruskalwallis}=N.A.';
    %             end
    %         end
    % 
    %         if iscell(OldTitle{idx})
    %             newtitle = OldTitle{idx};
    %         else
    %             newtitle = OldTitle(idx);
    %         end
    %         newtitle{2} = b;
    %         a.Title.String = newtitle;
    %     else
    %         a.Title.String = OldTitle{idx};
    %         OldTitle{idx} = [];
    %     end
    % end

     % function i_savedata_thistab(~, ~)
     %     [~,idx]=ismember(focalg, tabnamelist);
     %     thisy = y{idx};
     %     T = table(thisy(:), thisc(:));
     %     T.Properties.VariableNames = {'ScoreLevel', 'GroupID'};
     %     %T=sortrows(T,'ScoreLevel','descend');
     %     %T=sortrows(T,'GroupID');
     %     gui.i_exporttable(T, true, 'Tviolindata','ViolinPlotTable');
     % end


    function in_callback_savedata_alltab(~, ~)
%         [~,idx]=ismember(focalg, tabnamelist);
%         thisy = y{idx};
     
        T = table();
        for tabidx=1:n
            % g = tabnamelist(tabidx);
            thisy = y{tabidx};
            t = table(thisy(:));
            t.Properties.VariableNames = matlab.lang.makeValidName(tabnamelist(tabidx));
            T = [T, t];
        end
         t = table(thisc(:));
         t.Properties.VariableNames = {'GroupID'};
         T = [t, T];
         %T=sortrows(T,'ScoreLevel','descend');
         %T=sortrows(T,'GroupID');
         gui.i_exporttable(T, true, 'Tviolindata','ViolinPlotTable');
     end
 

    function in_callback_savedata(~, ~)
        [~, idxlabel]= findgroups(string(thisc(:)));
        T=table();
        for tabidx=1:n
            g = tabnamelist(tabidx);
            thisy = y{tabidx};
            
            a1=grpstats(thisy, thisc(:), @mean);
            a2=grpstats(thisy, thisc(:), @median);
            %{
            assignin("base","thisy",thisy);
            assignin("base","thisc",thisc);
            a1=splitapply(@mean, thisy, thisc(:));
            a2=splitapply(@median, thisy, thisc(:));
            %}
            
            t = table(a1, a2);
            t.Properties.RowNames = idxlabel;
            t.Properties.VariableNames = matlab.lang.makeValidName({sprintf('Mean_%s',g), sprintf('Median_%s',g)});
            T = [T, t];
        end
        T = rows2vars(T);
        gui.i_exporttable(T, true, 'Tviolindata','ViolinPlotTable');
    end

end
