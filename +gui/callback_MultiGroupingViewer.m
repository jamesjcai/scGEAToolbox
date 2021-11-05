function callback_MultiGroupingViewer(src,~)
    FigureHandle=src.Parent.Parent;
    sce=guidata(FigureHandle);    
    if isempty(sce.c_cluster_id) ||...
       size(sce.s,1)~=length(sce.c_cluster_id) ||...
       isempty(sce.c_batch_id) ||...
       size(sce.s,1)~=length(sce.c_batch_id)
       warndlg('This function requires both BATCH_ID and CLUSTER_ID are defined.');
       return;
    end
    gui.sc_multigroupings(sce,'Batch','Cluster');
end
