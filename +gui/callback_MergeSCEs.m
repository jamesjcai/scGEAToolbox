function [requirerefresh,s]=callback_MergeSCEs(src)
    requirerefresh=false;
    s='';
    a=evalin('base','whos');
    b=struct2cell(a);
    valididx=ismember(b(4,:),'SingleCellExperiment');
    if sum(valididx)<1
        warndlg('No SCE variables in Workspace.');
        return;        
    elseif sum(valididx)<2
        warndlg('Need at least two SCEs in Workspace.');
        return;
    end
    
    b=b(:,valididx);
    a=a(valididx);

    answer = questdlg('Current SCE will be replaced. Continue?');
    if ~strcmp(answer, 'Yes'), return; end
    FigureHandle=src.Parent.Parent;
    %sce=guidata(FigureHandle);
    
    [indx,tf]=listdlg('PromptString',{'Select SCEs:'},...
        'liststring',b(1,:),'SelectionMode','multiple');
    if tf==1
        if length(indx)<2
            warndlg('Need at least two selected SCEs.');
            return;
        end
        
        answer = questdlg('Which set operation method to merge genes?', 'Merging method',...
                          'Intersect', 'Union', 'Intersect');
        if ~ismember(answer, {'Union', 'Intersect'}), return; end
        methodtag = lower(answer);
        
        insce=cell(1,length(indx));
        s="";
        for k=1:length(indx)
            insce{k}=evalin('base',a(indx(k)).name);
            s=sprintf('%s,%s',s,a(indx(k)).name);
        end
        s=s(2:end);
        fprintf('>> sce=sc_mergesces({%s},''%s'');\n',s,methodtag);
        sce=sc_mergesces(insce,methodtag);
        guidata(FigureHandle,sce);
        requirerefresh=true;
    else
        return;
    end
end
