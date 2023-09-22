function xcallback_CrossTabulation_old(src, ~)
FigureHandle = src.Parent.Parent;
sce = guidata(FigureHandle);

[thisc1, clable1, thisc2, clable2] = gui.i_select2class(sce);
if isempty(thisc1) || isempty(thisc2)
    return;
end

%{
[thisc1,clable1]=gui.i_select1class(sce);
if isempty(thisc1), return; end
uiwait(helpdlg(sprintf('First grouping varible (%s) selected.',clable1)));
[thisc2,clable2]=gui.i_select1class(sce);
if isempty(thisc2), return; end
uiwait(helpdlg(sprintf('Second grouping varible (%s) selected.',clable2)));
% if strcmp(clable1,clable2), return; end
%}

answer = questdlg('Show groups by?', ...
    'Sorted Variable', ...
    clable1, clable2, clable1);
switch answer
    case clable1
        %[thisc1,thisc2]=i_sortc(thisc1,thisc2);
        thiscA = thisc1;
        thiscB = thisc2;
        clabel = clable1;
        llabel = clable2;
    case clable2
        %[thisc1,thisc2]=i_sortc(thisc2,thisc1);
        thiscA = thisc2;
        thiscB = thisc1;
        clabel = clable2;
        llabel = clable1;
        % case 'No sort'
    otherwise
        return;
end
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
%[~,cL]=grp2idx(thisc2);
%legend(cL);
% lgd=legend({'see below'},'Location','bestoutside');
% title(lgd,clable2)
% lgd=legend(labelsy,'Location','bestoutside');
% title(lgd,llabel);

subplot(212)
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
lgd = legend(labelsy, 'Location', 'bestoutside');
title(lgd, llabel);


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

    function i_sortbymean(~, ~)

        if sizesorted
            [~, idx] = sort(a, 'ascend');
            sizesorted = false;
        else
            [~, idx] = sort(a, 'descend');
            sizesorted = true;
        end
        cLx_sorted = cLx(idx);

        %[~,cL,noanswer]=gui.i_reordergroups(thisc,cLx_sorted);
        %if noanswer, return; end
        cla
        pkg.i_violinplot(y, thisc, colorit, cLx_sorted);
end

        function i_saveCrossTable(~, ~)
            gui.i_exporttable(T, false, 'TCrosstab');
            %         if ~(ismcc || isdeployed)
            %             labels = {'Save Cross-table to variable named:'};
            %             vars = {'TCrosstab'};
            %             values = {T};
            %             export2wsdlg(labels,vars,values);
            %         else
            %             gui.i_exporttable(T,false,'TCrosstab');
            %         end
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
