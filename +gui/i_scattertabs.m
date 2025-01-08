function i_scattertabs(y, tabnamelist, thisx, xlabelv, parentfig)
%see also: gui.i_violinplot, gui.sc_uitabgrpfig_vioplot

%assignin('base',"y",y)
%assignin('base',"thisx",thisx)

if nargin<4, parentfig = []; end
tabnamelist = string(tabnamelist);

%[~, cLorder]=grp2idx(thisx);
xlabelv = strrep(xlabelv, '_', '\_');

fw = gui.gui_waitbar;
% isdescend = false;

% thisx = strrep(string(thisx), '_', '\_');
% colorit = true;


import mlreportgen.ppt.*;
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, '..', 'resources', 'Misc', 'myTemplate.pptx');

hx = gui.myFigure;
hFig=hx.FigureHandle;
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
h0=cell(n,1);

OldTitle = cell(n,1);
for k = 1:n
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',tabnamelist(k)));
    ax0{k} = axes('parent',tab{k});
    h0{k} = scatter(ax0{k}, thisx(:), y{k}(:));
    xlabel(ax0{k}, xlabelv);
    ylabel(ax0{k}, strrep(tabnamelist(k), '_', '\_'));
    % pkg.i_violinplot(y{k}, thisx, true, cLorder);
    title(ax0{k}, strrep(tabnamelist(k), '_', '\_'));
end
  
tabgp.SelectionChangedFcn=@displaySelection;

hx.addCustomButton('off', @in_savedata, 'export.gif', 'Export data...');
hx.addCustomButton('off', @in_addregress, 'plotpicker-renko.gif', 'Add Regression Line...');
hx.addCustomButton('off', @in_addlocfitx, 'plotpicker-renkox.gif', 'Add Locfit Local Regression...');

hx.addCustomButton('off', @in_addlocfit, 'plotpicker-renko.gif', 'Add Locfit Local Regression...');


% hx.addCustomButton('off', @i_addsamplesize, "icon-mat-blur-linear-10.gif", 'Add Sample Size');
hx.addCustomButton('off', @i_savemainfig, "powerpoint.gif", 'Save Figure to PowerPoint File...');
hx.addCustomButton('off', @i_savemainfigx, "xpowerpoint.gif", 'Save Figure as Graphic File...');

%hx.addCustomButton('off', @i_invertcolor, "plotpicker-pie.gif", 'Switch BW/Color');
%hx.addCustomButton('off', @i_reordersamples, "plotpicker-errorbar.gif", 'Reorder Samples');

%hx.addCustomButton('off', @i_selectsamples, "plotpicker-errorbarx.gif", 'Select Samples');
%hx.addCustomButton('off', @i_sortbymean, "plotpicker-cra.gif", 'Sort Samples by Median');
hx.addCustomButton('on', @in_PickPlotMarker, 'plotpicker-rose.gif', 'Switch scatter plot marker type');
hx.addCustomButton('off', @in_BoxOnOff, 'RectGate.gif', 'Switch box on/off');

ybox = false;

if ~(ismcc || isdeployed)
    if license('test','Curve_Fitting_Toolbox') && ~isempty(which('curveFitter'))
        hx.addCustomButton('on', @in_curvefitter, 'icon-fa-stack-exchange-10.gif', 'Invoke curveFitter');
    end
end

gui.gui_waitbar(fw);
hx.show(parentfig);

     function in_curvefitter(~, ~)
         if ~(ismcc || isdeployed)
            [~,tabidx]=ismember(focalg, tabnamelist);            
            thisy = y{tabidx};
            curveFitter(thisx(:), thisy(:));
         end
     end
         
    function i_savemainfigx(~,~)
        [~,idx]=ismember(focalg, tabnamelist);     
        filter = {'*.jpg'; '*.png'; '*.tif'; '*.pdf'; '*.eps'};
        [filename, filepath] = uiputfile(filter,'Save Violin Plot', ...
            sprintf('ViolinPlot_%s',focalg));
        if ischar(filename)
            exportgraphics(ax0{idx}, [filepath, filename]);
        end
    end

    function i_savemainfig(~,~)
        answer = questdlg('Export to PowerPoint?');
        if ~strcmp(answer,'Yes'), return; end

        fw=gui.gui_waitbar_adv;
            OUTppt = [tempname, '.pptx'];
            ppt = Presentation(OUTppt, pth);
            open(ppt);
            images=cell(n,1);
            warning off
        for kx=1:n
            gui.gui_waitbar_adv(fw,kx./n,"Processing "+tabnamelist(kx)+" ...");
            images{kx} = [tempname, '.png'];
            tabgp.SelectedTab=tab{kx};
            saveas(tab{kx},images{kx});
            slide3 = add(ppt, 'Small Title and Content');
            replace(slide3, 'Title', tabnamelist(kx));
            replace(slide3, 'Content', Picture(images{kx}));        
        end
            close(ppt);
            rptview(ppt);      
            gui.gui_waitbar_adv(fw);
    end

    function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        % disp("Viewing gene " + txt);
        [~,idx]=ismember(txt,tabnamelist);
        focalg = tabnamelist(idx);
    end

    function in_BoxOnOff(~,~)
        for tabidx=1:n
            if ~ybox
                box(ax0{tabidx},'on');
            else
                box(ax0{tabidx},'off');
            end
            ybox = ~ybox;
        end        
    end

    function in_PickPlotMarker(~, ~)
        %answer = questdlg('Box on?','');
        %ybox = false;
        %if isempty(answer), return; end
        %if strcmp(answer,'Yes'), ybox = true; end
        s1 = 10 * randi(10);
        s2 = 50 * randi(10);        
        for tabidx=1:n
            if h0{tabidx}.Marker == '.'
                h0{tabidx}.Marker = 'o';
                h0{tabidx}.SizeData = s1;
            else
                h0{tabidx}.Marker = '.';
                h0{tabidx}.SizeData = s2;
            end
            if ybox
                box(ax0{tabidx},'on');
            else
                box(ax0{tabidx},'off');
            end
        end
    end

    function in_addlocfit(~, ~)
        for tabidx=1:n
            tabgp.SelectedTab=tab{tabidx};
            thisax = ax0{tabidx};
            thisy = y{tabidx};
            %[~,tabidx]=ismember(focalg, tabnamelist);
            %thisy = y{tabidx};
            [y_fit] = pkg.e_locfit(thisy(:), thisx(:));
            %[y_fit2] = malowess(thisx(:), thisy(:));

            if license('test','Curve_Fitting_Toolbox') && ~isempty(which('curveFitter'))
                y_fit2 = smooth(thisx(:), thisy(:), 0.5, 'lowess');
                %y_fit2 = fit([thisx(:), thisy(:)], 'lowess', 'Span', 0.2);
                % assignin("base",'y',thisy(:))
                % assignin("base",'x',thisx(:))
                % assignin("base",'y_fit1',y_fit)
                %assignin("base",'y_fit2',y_fit2)
            end
            

            %thisax = ax0{tabidx};
            hold(thisax,"on");
            [sortedx, idxx]=sort(thisx(:));
            plot(thisax, sortedx, y_fit(idxx), '-','LineWidth', 2);
            plot(thisax, sortedx, y_fit2(idxx), '-','LineWidth', 2);
            hold(thisax,"off");
        end
        [~,tabidx]=ismember(focalg, tabnamelist);
        tabgp.SelectedTab=tab{tabidx};
    end

    function in_addlocfitx(~, ~)
        % [~,tabidx]=ismember(focalg, tabnamelist);
        % thisy = y{tabidx};
        % [y_fit] = pkg.e_locfit(thisy(:), thisx(:));
        % thisax = ax0{tabidx};
        % hold(thisax,"on");
        % [sortedx, idxx]=sort(thisx(:));
        % plot(thisax, sortedx, y_fit(idxx), '-','LineWidth', 2);
        % hold(thisax,"off");

        [idxx] = gui.i_selmultidlg(tabnamelist, tabnamelist, hFig);
        if isempty(idxx), return; end
        if idxx == 0, return; end
        
        tabnamelist_a = tabnamelist(idxx);

        Y = [];
        for k=1:n
            if ismember(k, idxx)
                Y = [Y y{k}(:)];
            end
        end
        [Y_fit] = pkg.e_locfit(Y, thisx(:));
        answer = questdlg('ZScore transform data?','');
        switch answer
            case 'Yes'
                zs = true;
            case 'No'
                zs = false;
            otherwise
                return;
        end
        f = figure;
        f.Position(3)=f.Position(3)*1.8;
        hold on
        [sortedx, idxx]=sort(thisx(:));
        Y_fit = Y_fit(idxx,:);
        
        Pk = [];
        for k = 1:size(Y_fit, 2)
            if zs
                Pk = [Pk, plot(sortedx, zscore(Y_fit(:, k)), '-', 'LineWidth', 3)];
            else
                plot(thisx(:), Y(:,k), '.', 'markersize', 1);
                Pk = [Pk, plot(sortedx, Y_fit(:, k), '-', 'LineWidth', 3)];
            end
        end
        box on
        tabnamelist_a = strrep(tabnamelist_a,'_','\_');
        legend(Pk, tabnamelist_a, 'location', 'eastoutside');
        xlim([0, max(thisx(:))]);        
    end

    function in_addregress(~, ~)
        for tabidx=1:n
            tabgp.SelectedTab=tab{tabidx};
            thisax = ax0{tabidx};
            thisy = y{tabidx};
            %a = hFig.get("CurrentAxes");
            if isempty(OldTitle{tabidx})
                OldTitle{tabidx} = thisax.Title.String;
                if size(thisy, 2) ~= length(thisx)
                    thisy = thisy.';
                end
                % tbl = pkg.e_grptest(thisy, thisx);
                A  = fitlm(thisx(:), thisy(:));
                b = sprintf('Coeff. = %.2f; Adj. R^2 = %.2f; p-value = %.2e',...
                    A.Coefficients.Estimate(2), A.Rsquared.Adjusted, A.coefTest);
                % tbl = A.anova;
                % if ~isempty(tbl) && istable(tbl)
                %     b = sprintf('p_{anova} = %.2e', tbl.pValue);
                % else
                %         b='p_{anova} = N.A.';
                % end

                DM = [thisx(:), ones(size(thisx(:)))];
                B = DM \ thisy(:);
                y1 = DM * B;

                p = polyfit(thisx(:), thisy(:), 2);
                y_fit = polyval(p,thisx(:));

                hold(thisax,"on");
                [sortedx, idxx]=sort(thisx(:));
                plot(thisax, sortedx, y1(idxx), '-', 'LineWidth', 1);
                plot(thisax, sortedx, y_fit(idxx), '-','LineWidth', 1);
                hold(thisax,"off");

                if iscell(OldTitle{tabidx})
                    newtitle = OldTitle{tabidx};
                else
                    newtitle = OldTitle(tabidx);
                end
                newtitle{2} = b;
                thisax.Title.String = newtitle;
            else
                thisax.Title.String = OldTitle{tabidx};
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

    % function i_savedata(~, ~)
    %     [~,idx]=ismember(focalg, tabnamelist);
    %     thisy = y{idx};
    % 
    %     T = table(thisy(:), thisc(:));
    %     T.Properties.VariableNames = {'ScoreLevel', 'GroupID'};
    %     %T=sortrows(T,'ScoreLevel','descend');
    %     %T=sortrows(T,'GroupID');
    %     gui.i_exporttable(T, true, 'Tviolindata','ViolinPlotTable');
    % end

    function in_savedata(~, ~)        
        T=table(thisx);
        T.Properties.VariableNames = matlab.lang.makeValidName({xlabelv});
        for tabidx=1:n            
            thisy = y{tabidx};
            t = table(thisy(:));            
            t.Properties.VariableNames = matlab.lang.makeValidName(tabnamelist(tabidx));
            T = [T, t];
        end        
        gui.i_exporttable(T, true, 'Tscatterdata','ScatterPlotTable');
    end
end
