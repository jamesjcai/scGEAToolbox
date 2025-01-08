function callback_CrossTabulation(src, ~)

FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[thisc1, clabel1, thisc2, clabel2] = gui.i_select2class(sce, true);


% [c, cL, noanswer] = gui.i_reordergroups(thisc1, [], FigureHandle);
% thisc1=cL(c);
% [c, cL, noanswer] = gui.i_reordergroups(thisc2, [], FigureHandle);
% thisc2=cL(c);

if isempty(thisc1), return; end

if isempty(thisc2)
    hx = gui.myFigure;
    T=tabulate(thisc1);
    y = T(:,2);
    if iscell(y), y = cell2mat(y); end
    bar(y,'FaceColor', "flat");    
    colormap("turbo")
    labelsx=string(T(:,1));

    xticks(1:length(labelsx));
    labelsx1 = strrep(labelsx, '_', '\_');
    xticklabels(labelsx1);
    hx.show(FigureHandle);
    return;
end




fw = gui.gui_waitbar;

% sizesorted = false;
labelsx='';
labelsy='';
T=[];

hx = gui.myFigure;
tab=cell(2,1);
ax0=cell(2,1);
ax=cell(2,2);

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

hx.addCustomButton('off', @i_saveCrossTable, "export.gif", 'Save cross-table');
gui.gui_waitbar(fw);
hx.show(FigureHandle);


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
        % assignin("base","y",y);
        b = bar(y, 'stacked', 'FaceColor', "flat");
        %colormap(prism(size(y,2)));
        colormap(turbo);
        try
        for kx = 1:size(y, 2)
            b(kx).CData = kx;
        end
        catch
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
        try
        for kx = 1:size(y, 2)
            b(kx).CData = kx;
        end
        catch
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

    function i_saveCrossTable(~, ~)
        gui.i_exporttable(T, false, 'Tcrosstabul','CrosstabulTable');
    end
end
