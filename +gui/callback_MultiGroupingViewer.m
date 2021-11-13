function callback_MultiGroupingViewer(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);
        
    [thisc1,clable1,thisc2,clable2]=gui.i_select2class(sce);
     if isempty(thisc1) || isempty(thisc2)
         return;
     end
%     if isempty(sce.c_cluster_id) ||...
%        size(sce.s,1)~=length(sce.c_cluster_id) ||...
%        isempty(sce.c_batch_id) ||...
%        size(sce.s,1)~=length(sce.c_batch_id)
%        warndlg('This function requires both BATCH_ID and CLUSTER_ID are defined.');
%        return;
%     end
    
    gui.sc_multigroupings(sce,thisc1,thisc2,clable1,clable2);
end
