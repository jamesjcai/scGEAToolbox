function callback_CrossTabulation(src,~)
    FigureHandle=src.Parent.Parent;    
    sce=guidata(FigureHandle);  
    
    [thisc1,clable1,thisc2,clable2]=gui.i_select2class(sce);
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
     
answer = questdlg('Sort by?', ...
	'Sorted Variable', ...
	clable1,clable2,'No sort','No sort');
switch answer
    case clable1
        [thisc1,thisc2]=i_sortc(thisc1,thisc2);
    case clable2
        [thisc2,thisc1]=i_sortc(thisc2,thisc1);
    case 'No sort'
    otherwise
        return;
end
    
    [T,~,~,labelsxy]=crosstab(thisc1,thisc2);


labelsx=labelsxy(:,1);
labelsx=labelsx(~cellfun('isempty',labelsx));
labelsy=labelsxy(:,2);
labelsy=labelsy(~cellfun('isempty',labelsy));

figure;
subplot(211)
bar(T,'stacked')
xticks(1:length(labelsx));
xticklabels(labelsx);
xlabel(clable1)
ylabel('# of cells')
%[~,cL]=grp2idx(thisc2);
%legend(cL);
% lgd=legend({'see below'},'Location','bestoutside');
% title(lgd,clable2)

subplot(212)
bar(T./sum(T,2),'stacked')
xlabel(clable1)
ylabel('% of cells')
% title(clable2); 
xticks(1:length(labelsx));
xticklabels(labelsx);
ylim([0 1]);
lgd=legend(labelsy,'Location','bestoutside');
title(lgd,clable2);


    labels = {'Save Cross-table to variable named:'};
    vars = {'TCrosstab'};
    values = {T};
    export2wsdlg(labels,vars,values);
end

function [thisc1,thisc2]=i_sortc(thisc1,thisc2)
        [~,idx]=unique(thisc1);
        thisc1a=thisc1(idx);
        thisc1b=thisc1;
        thisc1b(idx)=[];
        thisc1=[thisc1a; thisc1b];

        thisc2a=thisc2(idx);
        thisc2b=thisc2;
        thisc2b(idx)=[];
        thisc2=[thisc2a; thisc2b];
end
