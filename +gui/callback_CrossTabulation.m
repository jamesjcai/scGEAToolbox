function callback_CrossTabulation(src,~)
    FigureHandle=src.Parent.Parent;    
    sce=guidata(FigureHandle);  
    [thisc1,c1txt]=gui.i_select1class(sce);
    if isempty(thisc1), return; end
    uiwait(helpdlg(sprintf('First grouping varible (%s) selected.',c1txt)));
    [thisc2,c2txt]=gui.i_select1class(sce);
    if isempty(thisc2), return; end
    uiwait(helpdlg(sprintf('Second grouping varible (%s) selected.',c2txt)));
    if strcmp(c1txt,c2txt), return; end
    [T,~,~,labelsxy]=crosstab(thisc1,thisc2);


labelsx=labelsxy(:,1);
labelsx=labelsx(~cellfun('isempty',labelsx));
labelsy=labelsxy(:,2);
labelsy=labelsy(~cellfun('isempty',labelsy));

figure;
subplot(211)
bar(T,'stacked')
xticklabels(labelsx);
xlabel(c1txt)
ylabel('# of cells')
%[~,cL]=grp2idx(thisc2);
%legend(cL);
lgd=legend(labelsy,'Location','bestoutside');
title(lgd,c2txt)
subplot(212)
bar(T./sum(T,2),'stacked')
xlabel(c1txt)
ylabel('% of cells')
% title(c2txt); 
xticklabels(labelsx);
lgd=legend(labelsy,'Location','bestoutside');
title(lgd,c2txt);

    labels = {'Save Cross-table to variable named:'};
    vars = {'TCrosstab'};
    values = {T};
    export2wsdlg(labels,vars,values);
end


