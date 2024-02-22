function sc_uitabgrpfig_violin(sce, glist, thisc, parentfig)


if nargin<4, parentfig = []; end


[~, cLorder]=grp2idx(thisc);
cLorder = strrep(cLorder, '_', '\_');
[Xt] = gui.i_transformx(sce.X);

isdescend = false;
OldTitle = [];
thisc = strrep(string(thisc), '_', '\_');
colorit = true;



if ~isempty(parentfig) && isa(parentfig,'matlab.ui.Figure') 
    p = parentfig.Position;
    cx = [p(1)+p(3)/2 p(2)+p(4)/2];
else
    cx = [];
end

import mlreportgen.ppt.*;
pw1 = fileparts(mfilename('fullpath'));
pth = fullfile(pw1, 'resources', 'myTemplate.pptx');

hFig = figure("Visible","off");
hFig.Position(3) = hFig.Position(3) * 1.8;

if ~isempty(cx)
    px = hFig.Position;
    px_new = [cx(1)-px(3)/2 cx(2)-px(4)/2];
else
    px_new = [];
end

n = length(glist);
% a = getpref('scgeatoolbox', 'prefcolormapname', 'autumn');

tabgp = uitabgroup();

idx = 1;
focalg = glist(idx);
tab=cell(n,1);
ax0=cell(n,1);

for k=1:n
    % c = sce.X(sce.g == glist(k), :);
    % if issparse(c), c = full(c); end
    tab{k} = uitab(tabgp, 'Title', sprintf('%s',glist(k)));
    ax0{k} = axes('parent',tab{k});
    y = full(Xt(upper(sce.g) == upper(glist(k)), :));
    pkg.i_violinplot(y, thisc, true, cLorder);
    title(ax0{k}, strrep(glist(k), '_', '\_'));
    % subtitle(ax0{k}, gui.i_getsubtitle(c));
    % gui.i_setautumncolor(c, a, true, any(c==0));
end
  

tabgp.SelectionChangedFcn=@displaySelection;
tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
uipushtool(tb, 'Separator', 'off');
pkg.i_addbutton2fig(tb, 'on',  @i_genecards, 'fvtool_fdalinkbutton.gif', 'GeneCards...');
% pkg.i_addbutton2fig(tb, 'on', {@i_PickColorMap, c}, 'plotpicker-compass.gif', 'Pick new color map...');

pkg.i_addbutton2fig(tb, 'off', @i_savedata, 'export.gif', 'Export data...');
pkg.i_addbutton2fig(tb, 'off', @i_testdata, 'icon-fa-stack-exchange-10.gif', 'ANOVA/T-test...');
pkg.i_addbutton2fig(tb, 'off', @i_addsamplesize, "icon-mat-blur-linear-10.gif", 'Add Sample Size');
pkg.i_addbutton2fig(tb, 'off', @i_savemainfig, "powerpoint.gif", 'Save Figure to PowerPoint File...');
pkg.i_addbutton2fig(tb, 'off', @i_invertcolor, "plotpicker-pie.gif", 'Switch BW/Color');
pkg.i_addbutton2fig(tb, 'off', @i_reordersamples, "plotpicker-errorbar.gif", 'Reorder Samples');

pkg.i_addbutton2fig(tb, 'off', @i_selectsamples, "plotpicker-errorbarx.gif", 'Select Samples');
pkg.i_addbutton2fig(tb, 'off', @i_sortbymean, "plotpicker-cra.gif", 'Sort Samples by Median');
pkg.i_addbutton2fig(tb, 'off', @i_renametitle, "icon-mat-touch-app-10.gif", 'Change Plot Title');
%pkg.i_addbutton2fig(tb, 'on', @i_viewgenenames, 'HDF_point.gif', 'Show Gene Names');
pkg.i_addbutton2fig(tb, 'on', @i_resizewin, 'HDF_pointx.gif', 'Resize Plot Window');

if isempty(px_new)
    movegui(hFig,'center');
else    
    movegui(hFig, px_new);
end
drawnow;
hFig.Visible=true;


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
            gui.gui_waitbar_adv(fw,kx./n,"Processing "+glist(kx)+" ...");
            images{kx} = [tempname, '.png'];
            tabgp.SelectedTab=tab{kx};
            saveas(tab{kx},images{kx});
            slide3 = add(ppt, 'Small Title and Content');
            replace(slide3, 'Title', glist(kx));
            replace(slide3, 'Content', Picture(images{kx}));        
        end
            close(ppt);
            rptview(ppt);      
            gui.gui_waitbar_adv(fw);
    end

    % function i_linksubplots(~,~)        
    %     hlink = linkprop([ax{idx,1},ax{idx,2}],{'CameraPosition','CameraUpVector'});
    % end

    function displaySelection(~,event)
        t = event.NewValue;
        txt = t.Title;
        % disp("Viewing gene " + txt);
        [~,idx]=ismember(txt,glist);
        focalg = glist(idx);
    end

    function i_genecards(~, ~)
        web(sprintf('https://www.genecards.org/cgi-bin/carddisp.pl?gene=%s', focalg),'-new');
    end


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

        %[~,cL,noanswer]=gui.i_reordergroups(thisc,cLx_sorted);
        %if noanswer, return; end
        b = hFig.get("CurrentAxes");
        cla(b);
        cLorder = cLx_sorted;
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_reordersamples(~, ~)
        [~, cLorder, noanswer] = gui.i_reordergroups(thisc);

        % cLorder
        if noanswer, return; end
        b = hFig.get("CurrentAxes");
        cla(b);
        pkg.i_violinplot(y, thisc, colorit, cLorder);
    end


    function i_selectsamples(~, ~)
        [~,cLorder]=grp2idx(thisc);
        [newidx] = gui.i_selmultidlg(cLorder, cLorder);
        if isempty(newidx), return; end
        picked=ismember(thisc,cLorder(newidx));
%        [~, cLorder, noanswer] = gui.i_reordergroups(thisc);
%        % cLorder
%        if noanswer, return; end
        
        cLorder=cLorder(ismember(cLorder,cLorder(newidx)));
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
                b='p\_ttest=N.A.; p\_wilcoxon=N.A.';
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


    function i_savedata(~, ~)    
        T = table(y(:), thisc(:));
        T.Properties.VariableNames = {'ScoreLevel', 'GroupID'};
        %T=sortrows(T,'ScoreLevel','descend');
        %T=sortrows(T,'GroupID');
        gui.i_exporttable(T, true, 'Tviolindata','ViolinPlotTable');

    end

end
