function callback_CrossTabulation(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[thisc1, clable1, thisc2, clable2] = gui.i_select2class(sce);


% [c, cL, noanswer] = gui.i_reordergroups(thisc1);
% thisc1=cL(c);
% [c, cL, noanswer] = gui.i_reordergroups(thisc2);
% thisc2=cL(c);

if isempty(thisc1) || isempty(thisc2)
    return;
end

sizesorted = false;
labelsx='';
labelsy='';
T=[];

answer = questdlg('Show groups by?', ...
    'Sorted Variable', ...
    clable1, clable2, 'Both', 'Both');
switch answer
    case clable1
        %[thisc1,thisc2]=i_sortc(thisc1,thisc2);
        thiscA = thisc1;
        thiscB = thisc2;
        clabel = clable1;
        llabel = clable2;
        in_crossplot(thiscA,thiscB);
        return;
    case clable2
        %[thisc1,thisc2]=i_sortc(thisc2,thisc1);
        thiscA = thisc2;
        thiscB = thisc1;
        clabel = clable2;
        llabel = clable1;
        % case 'No sort'  
        in_crossplot(thiscA,thiscB);
        return;
    case 'Both'
        thiscA = thisc1;
        thiscB = thisc2;
        clabel = clable1;
        llabel = clable2;
        in_crossplot(thiscA,thiscB);
        thiscA = thisc2;
        thiscB = thisc1;
        clabel = clable2;
        llabel = clable1;
        % case 'No sort'  
        in_crossplot(thiscA,thiscB);
        return;
        
    otherwise
        return;
end



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


    f0 = figure('Visible', false);
    sizesorted = false;
    subplot(211)
    in_plot1;
    subplot(212)
    in_plot2;
    tb = uitoolbar(f0);
    pkg.i_addbutton2fig(tb, 'off', @i_saveCrossTable, "export.gif", 'Save cross-table');
    pkg.i_addbutton2fig(tb, 'off', {@gui.i_savemainfig, 3}, "powerpoint.gif", 'Save Figure to PowerPoint File...');
    pkg.i_addbutton2fig(tb, 'on', @gui.i_pickcolormap, 'plotpicker-compass.gif', 'Pick new color map...');
    pkg.i_addbutton2fig(tb, 'on', @gui.i_invertcolor, 'plotpicker-comet.gif', 'Invert colors');
    % pkg.i_addbutton2fig(tb,'off',@i_reordersamples, ...
    %     "xpowerpoint.gif",'Reorder Samples');
    pkg.i_addbutton2fig(tb, 'off', @i_sortbymean, ...
        "xpowerpoint.gif", 'Sort Samples by Size');
    movegui(f0, 'center');
    set(f0, 'Visible', true);
    end


    function in_plot1
        y = T;
        b = bar(y, 'stacked', 'FaceColor', "flat");
        %colormap(prism(size(y,2)));
        colormap(turbo);
        for k = 1:size(y, 2)
            b(k).CData = k;
        end
        xticks(1:length(labelsx));
        xticklabels(labelsx);
        xlabel(clabel)
        ylabel('# of cells')
    end        

    function in_plot2
        y = T ./ sum(T, 2);
        b = bar(y, 'stacked', 'FaceColor', "flat");
        %colormap(prism(size(y,2)));
        colormap(turbo);
        for k = 1:size(y, 2)
            b(k).CData = k;
        end
        xlabel(clabel)
        ylabel('% of cells')
        % title(clable2);
        xticks(1:length(labelsx));
        xticklabels(labelsx);
        ylim([0, 1]);
        labelsy = strrep(labelsy, '_', '\_');
        lgd = legend(labelsy, 'Location', 'bestoutside');
        title(lgd, llabel);
    end        

    function i_sortbymean(~, ~)
        if ~sizesorted
            [~, idx] = sort(sum(T, 2), 'descend');
            T = T(idx, :);
            labelsx = labelsx(idx);
            sizesorted = true;
        else
            [idx] = gui.i_selmultidlg(labelsx, natsort(labelsx), FigureHandle);
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
        % "Tcellattrib","CellAttribTable"
        % "Tviolindata","ViolinPlotTable"
        % "Tcrosstabul","CrosstabulTable"
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
