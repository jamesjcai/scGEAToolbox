function [requirerefresh]=callback_MergeSCEs(src)
    requirerefresh=false;
    answer = questdlg('Current SCE will be replaced. Continue?');
    if ~strcmp(answer, 'Yes'), return; end    
    FigureHandle=src.Parent.Parent;
    %sce=guidata(FigureHandle);
    a=evalin('base','whos');
    b=struct2cell(a);
    valididx=ismember(b(4,:),'SingleCellExperiment');
    if sum(valididx)<2
        errordlg('Need at least two SCE variables in workspace.');
        return;
    end
    b=b(:,valididx);
    a=a(valididx);
    [indx,tf]=listdlg('PromptString',{'Select SCEs:'},...
        'liststring',b(1,:),'SelectionMode','multiple');
    if tf==1
        if length(indx)<2
            errordlg('Need at least two selected SCEs.');
            return;
        end
        insce=cell(1,length(indx));
        for k=1:length(indx)
            insce{k}=evalin('base',a(indx(k)).name);
        end
        sce=sc_mergesces(insce);
        guidata(FigureHandle,sce);
        requirerefresh=true;
    else
        return;
    end
end
