function callback_CrossTabulation(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[thisc1, clabel1, thisc2, clabel2] = gui.i_select2class(sce);


% [c, cL, noanswer] = gui.i_reordergroups(thisc1, [], FigureHandle);
% thisc1=cL(c);
% [c, cL, noanswer] = gui.i_reordergroups(thisc2, [], FigureHandle);
% thisc2=cL(c);

if isempty(thisc1) || isempty(thisc2)
    return;
end

fw = gui.gui_waitbar;

sizesorted = false;
labelsx='';
labelsy='';
T=[];

hFig=figure("Visible","off", "DockControls", "off");

[px_new] = gui.i_getchildpos(FigureHandle, hFig);


tabgp = uitabgroup();
for k=1:2
    switch k
        case 1
            thiscA = thisc1;
            thiscB = thisc2;
            clabel = clabel1;
            llabel = clabel2;
        case 2
            thiscA = thisc2;
            thiscB = thisc1;
            clabel = clabel2;
            llabel = clabel1;
    end
    in_crossplot(thiscA,thiscB);
    tab{k} = uitab(tabgp, 'Title', sprintf('Tab%d',k));
    %tab{k} = uitab(tabgp, 'Title', sprintf('%s-%s',clabel,llabel));
    ax0{k} = axes('parent',tab{k});
    ax{k,1} = subplot(2,1,1);    
    in_plot1;
    ax{k,2} = subplot(2,1,2);
    in_plot2;
end

tb = findall(hFig, 'Tag', 'FigureToolBar'); % get the figure's toolbar handle
    % tb = uitoolbar(hFig);
    uipushtool(tb, 'Separator', 'off');
    % pkg.i_addbutton2fig(tb, 'off', [], "IMG00107.GIF", " ");
    pkg.i_addbutton2fig(tb, 'off', @i_saveCrossTable, "export.gif", 'Save cross-table');
    pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
    pkg.i_addbutton2fig(tb, 'on', @gui.i_pickcolormap, 'plotpicker-compass.gif', 'Pick new color map...');
    pkg.i_addbutton2fig(tb, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
    % pkg.i_addbutton2fig(tb,'off',@i_reordersamples, ...
    %     "xpowerpoint.gif",'Reorder Samples');
    % pkg.i_addbutton2fig(tb, 'off', @i_sortbymean, ...
    %     "xpowerpoint.gif", 'Sort Samples by Size');
%    movegui(hFig, 'center');
%    set(hFig, 'Visible', true);

if any(px_new<0)
    movegui(hFig, 'center');
else
    movegui(hFig, px_new);
end
drawnow;
gui.gui_waitbar(fw);
hFig.Visible=true;



    function in_crossplot(thiscA,thiscB)

        t = table(thiscA, thiscB);
        t = sortrows(t, [1, 2]);
        thiscA = t.thiscA;
        thiscB = t.thiscB;
        
        [T, ~, ~, labelsxy] = crosstab(thiscA, thiscB);
        
        labelsx = labelsxy(:, 1);
        labelsx = labelsx(~cellfun('isempty', labelsx));
        labelsy = labelsxy(:, 2);
        labelsy = labelsy(~cellfun('isempty', labelsy));


    % f0 = figure('Visible', false);
    % sizesorted = false;
    % subplot(211)
    % in_plot1;
    % subplot(212)
    % in_plot2;
    
    end


    function in_plot1
        y = T;
        b = bar(y, 'stacked', 'FaceColor', "flat");
        %colormap(prism(size(y,2)));
        colormap(turbo);
        for kx = 1:size(y, 2)
            b(kx).CData = kx;
        end
        xticks(1:length(labelsx));
        labelsx1 = strrep(labelsx, '_', '\_');
        xticklabels(labelsx1);       
        
        xlabel(clabel)
        ylabel('# of cells')
        labelsy1 = strrep(labelsy, '_', '\_');
        lgd = legend(labelsy1, 'Location', 'bestoutside');
        title(lgd, llabel);

    end        

    function in_plot2
        y = T ./ sum(T, 2);
        b = bar(y, 'stacked', 'FaceColor', "flat");
        %colormap(prism(size(y,2)));
        colormap(turbo);
        for kx = 1:size(y, 2)
            b(kx).CData = kx;
        end
        xlabel(clabel)
        ylabel('% of cells')
        % title(clabel2);
        xticks(1:length(labelsx));

        labelsx2 = strrep(labelsx, '_', '\_');
        xticklabels(labelsx2);
        ylim([0, 1]);
        labelsy2 = strrep(labelsy, '_', '\_');
        lgd = legend(labelsy2, 'Location', 'bestoutside');
        title(lgd, llabel);
    end        

    function i_sortbymean(~, ~)
        if ~sizesorted
            [~, idx] = sort(sum(T, 2), 'descend');
            T = T(idx, :);
            labelsx = labelsx(idx);
            sizesorted = true;
        else
            [idx] = gui.i_selmultidlg(labelsx, natsort(labelsx));
            % [~, idx] = sort(labelsx);
            T = T(idx, :);
            labelsx = labelsx(idx);
            sizesorted = false;
        end
        subplot(211)
        cla
        in_plot1;
        subplot(212)
        cla
        in_plot2;
    end

    function i_saveCrossTable(~, ~)
        gui.i_exporttable(T, false, 'Tcrosstabul','CrosstabulTable');
    end


end






                % function [thisc1,thisc2]=i_sortc(thisc1,thisc2)
                %         [~,idx]=unique(thisc1);
                %         thisc1a=thisc1(idx);
                %         thisc1b=thisc1;
                %         thisc1b(idx)=[];
                %         thisc1=[thisc1a; thisc1b];
                %
                %         thisc2a=thisc2(idx);
                %         thisc2b=thisc2;
                %         thisc2b(idx)=[];
                %         thisc2=[thisc2a; thisc2b];
                % end
